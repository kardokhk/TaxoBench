"""
Evaluation of taxonomic classification results.

CSV-only output (no plots). No cross-validation logic.

Per-(level x metric) CSVs:
    Levels:  Class, Order, Family, Genus, Species
    Metrics: precision, recall, f-score, detection rate
"""

from os import stat, makedirs
from os.path import exists, join
from collections import defaultdict

import pandas as pd
from biom.parse import MetadataMap
from biom.exception import TableException
from biom import load_table
from biom.cli.util import write_biom_table
from biom.table import Table


# Mapping between human-readable level names and the integer level used by
# biom-style 7-level lineages (k/d, p, c, o, f, g, s)
LEVEL_NAME_TO_INT = {
    'Class':   2,
    'Order':   3,
    'Family':  4,
    'Genus':   5,
    'Species': 6,
}

DEFAULT_LEVELS = ['Class', 'Order', 'Family', 'Genus', 'Species']
DEFAULT_METRICS = ['precision', 'recall', 'f-score', 'detection_rate']


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_results(expected_table_fp, observed_table_fp,
                    feature_table_fp,
                    expected_taxa_fp, observed_taxa_fp,
                    taxonomy_level,
                    metadata,
                    out_file,
                    min_count,
                    sample_weight,
                    taxa_to_keep=None, md_key='taxonomy',
                    per_seq_precision=False, exclude=None):
    """Compute precision/recall/f-score and taxon-detection metrics for a
    single (mock x region x classifier x param) combination at one
    `taxonomy_level`.  Writes a long-format TSV to `out_file`.
    """
    results = []

    expected_table = mount_observations(
        expected_table_fp, min_count=min_count,
        taxonomy_level=taxonomy_level,
        taxa_to_keep=taxa_to_keep, md_key=md_key,
        normalize=False, clean_obs_ids=False, filter_obs=False)

    observed_table = mount_observations(
        observed_table_fp, min_count=min_count,
        taxonomy_level=taxonomy_level,
        taxa_to_keep=taxa_to_keep, md_key=md_key,
        normalize=False, clean_obs_ids=False, filter_obs=False)

    try:
        feature_table = load_table(feature_table_fp)
    except ValueError:
        raise ValueError("Couldn't parse BIOM table: {0}".format(feature_table_fp))

    for sample_id in feature_table.ids(axis="sample"):
        observed_sample_id = sample_id
        expected_sample_id = observed_sample_id.rsplit('-', 1)[0]
        try:
            n_exp_taxa, n_obs_taxa, tp, fp, fn, accuracy, detection = \
                compute_taxon_accuracy(observed_table, expected_table,
                                       actual_sample_id=observed_sample_id,
                                       expected_sample_id=expected_sample_id)
        except ZeroDivisionError:
            n_exp_taxa = n_obs_taxa = tp = fp = fn = accuracy = detection = -1.

        if per_seq_precision:
            n_exp_seq, n_obs_seq, precision, recall, fscore, per_seq_id, per_seq_match = \
                per_sequence_precision(
                    expected_taxa_fp, observed_taxa_fp, feature_table,
                    sample_id, taxonomy_level,
                    sample_weight=sample_weight, exclude=exclude)
        else:
            n_exp_seq = n_obs_seq = precision = recall = fscore = -1.
            per_seq_id = per_seq_match = "None"

        mock_id, region_id, method_id, param_id, taxa_level = metadata

        results.append((
            mock_id, region_id, method_id, param_id, taxa_level, sample_id,
            n_exp_taxa, n_obs_taxa, tp, fp, fn, accuracy, detection,
            n_exp_seq, n_obs_seq, precision, recall, fscore,
            per_seq_id, per_seq_match
        ))

    result = pd.DataFrame(results, columns=[
        "mock_name", "v_region", "method", "param", "taxa_level", "sample_id",
        "n_exp_taxa", "n_obs_taxa", "tp", "fp", "fn",
        "Taxon Accuracy Rate", "Taxon Detection Rate",
        "n_exp_seq", "n_obs_seq", "Precision", "Recall", "F-measure",
        "per_seq_ids", "per_seq_matches",
    ])
    result.to_csv(out_file, sep="\t", index=False, header=True)
    return result


def mount_observations(table_fp, min_count=0, taxonomy_level=6,
                       taxa_to_keep=None, md_key='taxonomy', normalize=True,
                       clean_obs_ids=True, filter_obs=True):
    """Load biom table, optionally filter, then collapse on taxonomy."""
    try:
        table = load_table(table_fp)
    except ValueError:
        raise ValueError("Couldn't parse BIOM table: {0}".format(table_fp))

    if filter_obs is True and min_count > 0 and taxa_to_keep is not None:
        try:
            table = filter_table(table, min_count, taxonomy_level,
                                 taxa_to_keep, md_key=md_key)
        except TableException:
            pass
        except TypeError:
            print("Missing taxonomic information in table " + table_fp)

        if table.is_empty():
            raise ValueError("Table is empty after filtering at {0}".format(table_fp))

    collapse_taxonomy = get_taxonomy_collapser(taxonomy_level, md_key=md_key)
    try:
        table = table.collapse(collapse_taxonomy, axis='observation',
                               min_group_size=1)
    except TableException:
        raise TableException("Failure to collapse taxonomy for table at: {0}".format(table_fp))
    except TypeError:
        raise TypeError("Failure to collapse taxonomy: {0}".format(table_fp))

    if normalize is True:
        table.norm(axis='sample')
    return table


def filter_table(table, min_count=0, taxonomy_level=None,
                 taxa_to_keep=None, md_key='taxonomy'):
    try:
        _taxa_to_keep = ';'.join(taxa_to_keep)
    except TypeError:
        _taxa_to_keep = None

    def f(data_vector, id_, metadata):
        taxonomy = metadata.get(md_key) if metadata else None
        if isinstance(taxonomy, str):
            taxonomy = [level.strip() for level in taxonomy.split(';')]
        elif isinstance(taxonomy, list):
            taxonomy = [level.strip() for level in taxonomy]
        else:
            taxonomy = None
        enough_levels = taxonomy_level is None or \
            (taxonomy is not None and len(taxonomy) >= taxonomy_level + 1)
        allowed_taxa = _taxa_to_keep is None or id_.startswith(_taxa_to_keep) or \
            (taxonomy is not None and ';'.join(taxonomy).startswith(_taxa_to_keep))
        sufficient_count = data_vector.sum() >= min_count
        return sufficient_count and allowed_taxa and enough_levels

    return table.filter(f, axis='observation', inplace=False)


def get_taxonomy_collapser(level, md_key='taxonomy', unassignable_label='unassigned'):
    def f(id_, md):
        try:
            taxonomy = md[md_key]
            if isinstance(taxonomy, str):
                levels = [l.strip() for l in taxonomy.split(';')]
            elif isinstance(taxonomy, list):
                if len(taxonomy) == 1 and ';' in taxonomy[0]:
                    levels = [l.strip() for l in taxonomy[0].split(';')]
                else:
                    levels = [l.strip() for l in taxonomy]
            else:
                levels = [unassignable_label]
        except (TypeError, KeyError):
            levels = [unassignable_label]

        while len(levels) <= level:
            levels.append('')

        clean_levels = []
        for i, level_str in enumerate(levels):
            prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
            if i < len(prefixes) and level_str.startswith(prefixes[i]):
                clean_levels.append(level_str[3:].strip())
            else:
                clean_levels.append(level_str.strip())

        for check_level in range(level, -1, -1):
            if check_level < len(clean_levels) and clean_levels[check_level]:
                if level == 5:
                    if check_level == 6 and clean_levels[6]:
                        species_parts = clean_levels[6].split()
                        return species_parts[0] if species_parts else unassignable_label
                    elif check_level == 5 and clean_levels[5]:
                        return clean_levels[5]
                    else:
                        return unassignable_label
                elif level == 6:
                    if check_level == 6 and clean_levels[6]:
                        return clean_levels[6]
                    elif check_level == 5 and clean_levels[5]:
                        return clean_levels[5]
                    else:
                        return unassignable_label
                else:
                    if check_level == level:
                        return clean_levels[level]
                    else:
                        return unassignable_label
        return unassignable_label
    return f


def compute_taxon_accuracy(actual_table, expected_table, actual_sample_id=None,
                           expected_sample_id=None, ws_strip=True):
    actual_obs_ids = get_observed_observation_ids(actual_table, actual_sample_id, ws_strip=ws_strip)
    expected_obs_ids = get_observed_observation_ids(expected_table, expected_sample_id, ws_strip=ws_strip)
    total_exp = len(expected_obs_ids)
    total_obs = len(actual_obs_ids)
    tp = len(actual_obs_ids & expected_obs_ids)
    fp = len(actual_obs_ids - expected_obs_ids)
    fn = len(expected_obs_ids - actual_obs_ids)
    if tp > 0:
        p = round((tp / (tp + fp)), 3)
        r = round((tp / (tp + fn)), 3)
    else:
        p, r = 0, 0
    return total_exp, total_obs, tp, fp, fn, p, r


def get_observed_observation_ids(table, sample_id=None, ws_strip=False):
    if sample_id is None:
        sample_id = table.ids(axis="sample")[0]
    result = []
    for observation_id in table.ids(axis="observation"):
        if table.get_value_by_ids(observation_id, sample_id) > 0.0:
            if ws_strip is True:
                observation_id = "".join(observation_id.split())
            result.append(observation_id)
    return set(result)


# ---------------------------------------------------------------------------
# Per-sequence precision / recall / f-score
# ---------------------------------------------------------------------------

def per_sequence_precision(expected_table_fp, actual_table_fp, feature_table,
                           sample_id, taxonomy_level, sample_weight=None, exclude=None):
    exp_fp = expected_table_fp
    if exists(exp_fp):
        if exists(actual_table_fp):
            obs_fp = actual_table_fp
        else:
            raise RuntimeError('taxonomy assignments do not exist for dataset {0}'.format(actual_table_fp))
        exp = observations_to_list(exp_fp, feature_table, sample_id)
        obs = observations_to_list(obs_fp, feature_table, sample_id)
        try:
            exp, obs, exp_ids = _multiple_match_kludge(exp, obs, fill_empty_observations=True)
        except AssertionError:
            print('AssertionError in:')
            print(actual_table_fp)
            raise
        exp_taxa, obs_taxa = load_prf(obs, exp, level=slice(0, taxonomy_level + 1), sort=False)
        if sample_weight:
            sample_weight = [feature_table.get_value_by_ids(line.split('\t')[0], sample_id) for line in exp]
        total_exp, total_obs, ps, rs, fs, detailed_output = compute_prf(
            exp_taxa, obs_taxa, sample_weight=sample_weight, exclude=exclude)
    else:
        total_exp, total_obs, ps, rs, fs, exp_ids, detailed_output = -1., -1., -1., -1., -1., "None", "None"
    return total_exp, total_obs, ps, rs, fs, exp_ids, detailed_output


def observations_to_list(obs_fp, actual_table, sample_id):
    obs = import_to_list(input_data=obs_fp, field=slice(None, 2), delim="\t")
    obs = [line for line in obs if
           actual_table.exists(line.split('\t')[0], "observation") and
           actual_table.get_value_by_ids(line.split('\t')[0], sample_id) != 0]
    return obs


def import_to_list(input_data, field=None, delim='\t'):
    if isinstance(input_data, str):
        with open(input_data, "r") as f:
            lines = [line.strip() for line in f]
    elif isinstance(input_data, list):
        lines = [line.strip() for line in input_data]
    else:
        raise ValueError("input_data must be a file path or a list of lines")

    if field is not None:
        result = [
            delim.join(line.split(delim)[field]) if isinstance(field, slice)
            else line.split(delim)[field]
            for line in lines
        ]
    else:
        result = lines
    return result


def _multiple_match_kludge(exp, obs, fill_empty_observations=True):
    obs = {i: t for i, t in [r.split('\t', 1) for r in obs]}
    exp_grouped = defaultdict(list)
    for exp_id, exp_taxon in [r.split('\t') for r in exp]:
        exp_grouped[exp_id].append(exp_taxon)
    if fill_empty_observations:
        for k in exp_grouped.keys():
            if k not in obs.keys():
                obs[k] = 'Unassigned'
    else:
        assert obs.keys() == exp_grouped.keys(), \
            'observed and expected read labels differ:\n' + \
            str(list(obs.keys())) + '\n' + str(list(exp_grouped.keys()))
    new_exp = []
    new_obs = []
    exp_ids = []
    for exp_id, exp_taxons in exp_grouped.items():
        obs_row = '\t'.join([exp_id, obs[exp_id]])
        for exp_taxon in exp_taxons:
            if obs[exp_id] == exp_taxon:
                row = '\t'.join([exp_id, exp_taxon])
                break
        else:
            row = '\t'.join([exp_id, exp_taxons[0]])
        new_exp.append(row)
        new_obs.append(obs_row)
        exp_ids.append(exp_id)
    return new_exp, new_obs, ",".join(exp_ids)


def load_prf(obs_fp, exp_fp, level=slice(0, 7), sort=False):
    exp_taxa = load_taxa(exp_fp, level=level, sort=sort)
    obs_taxa = load_taxa(obs_fp, level=level, sort=sort)
    if len(obs_taxa) != len(exp_taxa):
        raise RuntimeError('Lengths of expected and observed taxa do not match.')
    return exp_taxa, obs_taxa


def load_taxa(obs_fp, level=slice(0, 7), field=1, sort=False):
    raw_obs = import_to_list(obs_fp)
    if sort:
        raw_obs = sorted(raw_obs)
    obs = extract_taxa_names(raw_obs, field=field, level=level)
    obs = [';'.join([l.strip() for l in line.split(';')]) for line in obs
           if not line.startswith(('taxonomy', 'Taxon'))]
    return obs


def extract_taxa_names(infile, level=slice(6, 7), l_delim=';', field=None,
                       f_delim='\t', stripchars=None):
    line_list = import_to_list(infile, field=field, delim=f_delim)
    name_list = [l_delim.join(line.split(l_delim)[level]).strip(stripchars)
                 for line in line_list]
    return name_list


def compute_prf(exp, obs, sample_weight=None, exclude=None):
    """Compute precision, recall, F-measure for a 'mock' style comparison
    (one expected taxonomy per sequence; no cross-validation)."""
    total_exp, total_obs, p, r, f, detailed = precision_recall_fscore(
        exp, obs, sample_weight=sample_weight, exclude=exclude)
    return total_exp, total_obs, p, r, f, detailed


def precision_recall_fscore(exp, obs, sample_weight=None, exclude=None):
    if exclude is None:
        exclude = []
    if sample_weight is None or sample_weight is False:
        sample_weight = [1] * len(exp)
    total_exp = len(exp)
    total_obs = len(obs)
    detailed = []
    tp, fp, fn = 0, 0, 0
    for e, o, w in zip(exp, obs, sample_weight):
        if e in exclude:
            continue
        classification = evaluate_classification(o, e)
        detailed.append(classification)
        if classification == 'match':
            tp += w
        elif classification == 'underclassification':
            fn += w
        else:
            fp += w
            fn += w
    if tp == 0:
        return total_exp, total_obs, 0, 0, 0, ",".join(detailed)
    p = round((tp / (tp + fp)), 3)
    r = round((tp / (tp + fn)), 3)
    f = round((2 * p * r / (p + r)), 3)
    return total_exp, total_obs, p, r, f, ",".join(detailed)


def _expand_group_label(label, prefix=None):
    """Return the set of taxon names encoded by a (possibly grouped) label.

    Taxonova may emit a multi-value group label of the form
    ``[a|b|c]`` at ANY taxonomic rank when several taxa share an identical
    sub-region sequence (comment 2.9 fix, generalized to all ranks). For a
    plain single-taxon label we return a 1-element set; for a group label we
    return the full member set; for empty / ambiguous we return an empty
    set. Optional ``prefix`` (e.g. 's__', 'g__', 'f__') is stripped first.
    """
    if not label:
        return set()
    s = label.strip()
    if prefix and s.startswith(prefix):
        s = s[len(prefix):]
    if not s:
        return set()
    if s.startswith("[") and s.endswith("]") and "|" in s:
        return {p.strip() for p in s[1:-1].split("|") if p.strip()}
    return {s}


def _expand_species_group(label):
    """Backward-compatible alias for the species-rank case."""
    return _expand_group_label(label, prefix="s__")


_RANK_PREFIXES = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')


def _strip_rank_prefix(part):
    """Return (prefix, content) for a single taxon-string part."""
    p = part.strip()
    if p[:3] in _RANK_PREFIXES:
        return p[:3], p[3:]
    return '', p


def _ambiguous(s):
    return (not s) or s in ('AMBIGUOUS', 'Unclassified', 'Unassigned', 'No blast hit')


def evaluate_classification(obs_taxon, exp_taxon):
    """Score a (predicted, expected) pair as match / over /
    underclassification / misclassification.

    Taxonova multi-value group labels (``[a|b|c]``) at ANY rank are matched
    by set membership: a prediction counts as a match at rank L iff the
    intersection of the predicted and expected member sets at rank L is
    non-empty (comment-2.9 multi-rank generalization).
    """
    def normalize_unclassified(taxon_string):
        if not taxon_string:
            return taxon_string
        parts = [part.strip() for part in taxon_string.split(';')]
        normalized = []
        for part in parts:
            prefix, content = _strip_rank_prefix(part)
            if prefix:
                if content.startswith('unclassified_'):
                    normalized.append(prefix)
                else:
                    normalized.append(part)
            else:
                normalized.append(part)
        return '; '.join(normalized)

    normalized_obs = normalize_unclassified(obs_taxon)
    normalized_exp = normalize_unclassified(exp_taxon)

    # Fast path: literal equality (covers identical single-value strings).
    if normalized_obs == normalized_exp:
        return 'match'
    # Whole-string under/over heuristics still help for pure single-value
    # strings (group labels typically don't satisfy these).
    if normalized_obs in ('Unclassified', 'Unassigned', 'No blast hit'):
        return 'underclassification'

    obs_parts = [part.strip() for part in normalized_obs.split(';')]
    exp_parts = [part.strip() for part in normalized_exp.split(';')]

    # Pad parts to length 7 so missing ranks parse as empty.
    while len(obs_parts) < 7:
        obs_parts.append('')
    while len(exp_parts) < 7:
        exp_parts.append('')

    obs_sets = []
    exp_sets = []
    for i, prefix in enumerate(_RANK_PREFIXES):
        _, obs_content = _strip_rank_prefix(obs_parts[i])
        _, exp_content = _strip_rank_prefix(exp_parts[i])
        obs_sets.append(set() if _ambiguous(obs_content) else _expand_group_label(obs_content))
        exp_sets.append(set() if _ambiguous(exp_content) else _expand_group_label(exp_content))

    # Walk ranks from deepest (species) up to find the deepest expected level.
    # Then classify based on whether obs intersects exp at every shallower rank
    # AND at the deepest expected level.
    deepest_exp = -1
    for i in range(6, -1, -1):
        if exp_sets[i]:
            deepest_exp = i
            break
    if deepest_exp < 0:
        # expected has no rank info — treat as match
        return 'match'

    deepest_obs = -1
    for i in range(6, -1, -1):
        if obs_sets[i]:
            deepest_obs = i
            break

    if deepest_obs < 0:
        return 'underclassification'

    # Walk every rank up to min(deepest_exp, deepest_obs); if any rank
    # disagrees (non-empty exp at that rank with no intersection with obs)
    # → misclassification.
    for i in range(min(deepest_exp, deepest_obs) + 1):
        if exp_sets[i] and obs_sets[i]:
            if not (exp_sets[i] & obs_sets[i]):
                return 'misclassification'

    if deepest_obs < deepest_exp:
        return 'underclassification'
    if deepest_obs > deepest_exp:
        return 'overclassification'
    # Same depth and all rank intersections non-empty.
    return 'match'


# ---------------------------------------------------------------------------
# BIOM helpers
# ---------------------------------------------------------------------------

def add_metadata_to_biom_table(input_fp, taxonomy_map_fp, biom_output_fp, mode='observed'):
    if mode == 'observed':
        newbiom = load_table(input_fp)
        if stat(taxonomy_map_fp).st_size == 0:
            metadata = {}
        else:
            metadata = MetadataMap.from_file(taxonomy_map_fp,
                                             header=['Feature ID', 'taxonomy', 'c'])
        newbiom.add_metadata(metadata, 'observation')
        write_biom_table(newbiom, 'json', biom_output_fp)
    elif mode == 'expected':
        df = pd.read_csv(input_fp, sep="\t")
        df.set_index("Taxon", inplace=True)
        newbiom = Table(
            df.values,
            observation_ids=df.index.tolist(),
            sample_ids=df.columns.tolist(),
            observation_metadata=[{"taxonomy": idx} for idx in df.index]
        )
        write_biom_table(newbiom, 'json', biom_output_fp)


# ---------------------------------------------------------------------------
# Per-(level x metric) CSV writer
# ---------------------------------------------------------------------------

_METRIC_TO_COLUMN = {
    'precision':       'Precision',
    'recall':          'Recall',
    'f-score':         'F-measure',
    'fscore':          'F-measure',
    'detection_rate':  'Taxon Detection Rate',
}


def write_per_level_metric_csvs(long_df, out_dir,
                                levels=None, metrics=None,
                                row_key='method', col_key='param'):
    """Given a long-format dataframe with columns:
        mock_name, v_region, method, param, taxa_level, sample_id,
        Precision, Recall, F-measure, Taxon Detection Rate, ...

    write one CSV per (level x metric) into `out_dir`.

    File naming:  {out_dir}/{level}_{metric}.csv
    Each CSV is a wide table:  rows = `row_key` (e.g. method),
    columns = `col_key` (e.g. param), values = mean of the metric.
    """
    levels = levels or DEFAULT_LEVELS
    metrics = metrics or DEFAULT_METRICS
    makedirs(out_dir, exist_ok=True)

    written = []
    for lvl in levels:
        lvl_int = LEVEL_NAME_TO_INT.get(lvl)
        if lvl_int is None:
            continue
        sub = long_df[long_df['taxa_level'] == lvl_int]
        if sub.empty:
            sub = long_df[long_df['taxa_level'] == lvl]
        if sub.empty:
            continue
        for metric in metrics:
            col = _METRIC_TO_COLUMN.get(metric)
            if col is None or col not in sub.columns:
                continue
            try:
                pivot = sub.pivot_table(
                    index=row_key, columns=col_key, values=col, aggfunc='mean')
            except Exception:
                pivot = sub.groupby(row_key)[col].mean().to_frame()
            out_fp = join(out_dir, f"{lvl}_{metric}.csv")
            pivot.to_csv(out_fp)
            written.append(out_fp)
    return written
