from os import stat
from os.path import exists, join
import pandas as pd
from collections import defaultdict
from biom.parse import MetadataMap
from biom.exception import TableException
from biom import load_table
from biom.cli.util import write_biom_table
from biom.table import Table

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns



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
    """ Compute accuracy, detection, precision, recall, and f-measure 
    for result_tables at taxonomy_level

        taxonomy_level: taxonomic level to compute results
        per_seq_precision: bool
            Compute per-sequence precision/recall scores from expected
            taxonomy assignments?
        exclude: list
            taxonomies to explicitly exclude from precision scoring.
    """
    results = []
    # parse the expected table and collapse it at the specified taxonomic level
    expected_table = mount_observations(
        expected_table_fp,
        min_count=min_count,
        taxonomy_level=taxonomy_level,
        taxa_to_keep=taxa_to_keep,
        md_key=md_key, 
        normalize=False,
        clean_obs_ids=False,
        filter_obs=False)

    # parse the observed table and collapse it at the specified taxonomic level
    observed_table = mount_observations(
        observed_table_fp,
        min_count=min_count,
        taxonomy_level=taxonomy_level,
        taxa_to_keep=taxa_to_keep,
        md_key=md_key, 
        normalize=False,
        clean_obs_ids=False,
        filter_obs=False)
        

    # load the feature table without taxonomy assignment
    try:
        feature_table = load_table(feature_table_fp)
    except ValueError:
        raise ValueError("Couldn't parse BIOM table: {0}".format(feature_table_fp))

    # for each sample id
    for sample_id in feature_table.ids(axis="sample"):
        # when replicates simulated with ART, file name can be
        # expected = sample00
        # observed replicate(s) = sample00-1, sample00-2
        observed_sample_id = sample_id
        expected_sample_id = observed_sample_id.rsplit('-', 1)[0]
        # compute precision, recall, and f-measure
        try:
            n_exp_taxa, n_obs_taxa, tp, fp, fn, accuracy, detection = \
            compute_taxon_accuracy(
                observed_table, expected_table,
                actual_sample_id=observed_sample_id,
                expected_sample_id=expected_sample_id)
        except ZeroDivisionError:
            n_exp_taxa, n_obs_taxa, tp, fp, fn, accuracy, detection = \
                -1., -1., -1., -1., -1., -1., -1.

        # compute per-sequence precion / recall
        if per_seq_precision:
            n_exp_seq, n_obs_seq, precision, recall, fscore, per_seq_id, per_seq_match = \
                per_sequence_precision(
                    expected_taxa_fp, observed_taxa_fp, feature_table,
                    sample_id, taxonomy_level, sample_weight=sample_weight, exclude=exclude)
        else:
            n_exp_seq, n_obs_seq, precision, recall, fscore, per_seq_id, per_seq_match = \
                -1., -1., -1., -1., -1., "None", "None"
        
        # unpack metadata
        mock_id, region_id, method_id, param_id, taxa_level = metadata

        # log results
        results.append((
            mock_id, region_id, method_id, param_id, taxa_level, sample_id,
            n_exp_taxa, n_obs_taxa, tp, fp, fn, accuracy, detection,
            n_exp_seq, n_obs_seq, precision, recall, fscore, per_seq_id, per_seq_match
        ))

    # to dataframe
    result = pd.DataFrame(results, columns=["mock_name", "v_region", "method", "param", "taxa_level", "sample_id",
                                            "n_exp_taxa", "n_obs_taxa", "tp", "fp", "fn",
                                            "Taxon Accuracy Rate", "Taxon Detection Rate",
                                            "n_exp_seq", "n_obs_seq", "Precision", "Recall", "F-measure",
                                            "per_seq_ids", "per_seq_matches"
                                            ])
    
    # save to file
    result.to_csv(out_file, sep="\t", index=False, header=True)
    
    return result


def mount_observations(table_fp, min_count=0, taxonomy_level=6,
                       taxa_to_keep=None, md_key='taxonomy', normalize=True,
                       clean_obs_ids=True, filter_obs=True):
    '''load biom table, filter by abundance, collapse taxonomy, return biom.

    table_fp: path
        Input biom table.
    min_count: int
        Minimum abundance threshold; features detected at lower abundance are
        removed from table.
    taxonomy_level: int
        Taxonomic level at which to collapse table.
    taxa_to_keep: list of taxonomies to retain, others are removed before
        evaluation.
    md_key: str
        biom observation metadata key on which to collapse and filter.
    normalize: bool
        Normalize table to relative abundance across sample rows?
    clean_obs_ids: bool
        Remove '[]()' characters from observation ids? (these are removed from
        the ref db during filtering/cleaning steps, and should be removed from
        expected taxonomy files to avoid mismatches).
    filter_obs: bool
        Filter observations? filter_table will remove observations if taxonomy
        strings are shorter than taxonomic_level, count is less than min_count,
        or observation is not included in taxa_to_keep.
    '''
    try:
        table = load_table(table_fp)
    except ValueError:
        raise ValueError("Couldn't parse BIOM table: {0}".format(table_fp))

    if filter_obs is True and min_count > 0 and taxa_to_keep is not None:
        try:
            table = filter_table(table, min_count, taxonomy_level,
                                 taxa_to_keep, md_key=md_key)
        except TableException:
            # if all data is filtered out, move on to the next table
            pass

        except TypeError:
            print("Missing taxonomic information in table " + table_fp)

        if table.is_empty():
            raise ValueError("Table is empty after filtering at"
                             " {0}".format(table_fp))

    collapse_taxonomy = get_taxonomy_collapser(taxonomy_level, md_key=md_key)

    try:
        table = table.collapse(collapse_taxonomy, axis='observation',
                               min_group_size=1)
    except TableException:
        raise TableException("Failure to collapse taxonomy for table at:"
                             " {0}".format(table_fp))
    except TypeError:
        raise TypeError("Failure to collapse taxonomy: {0}".format(table_fp))

    if normalize is True:
        table.norm(axis='sample')

    return table

def filter_table(table, min_count=0, taxonomy_level=None,
                 taxa_to_keep=None, md_key='taxonomy'):
    # dynamically identiy and split if taxonomy is list instance or not
    # if filtering based on number of taxonomy levels, and this
    # observation has taxonomic information, and
    # there are a sufficient number of taxonomic levels
    # Table filtering here removed taxa that have insufficient levels
    try:
        _taxa_to_keep = ';'.join(taxa_to_keep)
    except TypeError:
        _taxa_to_keep = None

    def f(data_vector, id_, metadata):
        taxonomy = metadata.get(md_key) if metadata else None

        # Normalize taxonomy to a list if it's a string
        if isinstance(taxonomy, str):
            taxonomy = [level.strip() for level in taxonomy.split(';')]
        elif isinstance(taxonomy, list):
            taxonomy = [level.strip() for level in taxonomy]
        else:
            taxonomy = None

        # Check taxonomy level
        enough_levels = taxonomy_level is None or \
                        (taxonomy is not None and
                         len(taxonomy) >= taxonomy_level + 1)

        # Filter by allowed taxa
        allowed_taxa = _taxa_to_keep is None or \
            id_.startswith(_taxa_to_keep) or \
            (taxonomy is not None and
             ';'.join(taxonomy).startswith(_taxa_to_keep))

        # Filter by minimum count
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

        # Special handling for genus and species level comparisons
        # to handle SPINGO-style taxonomies with empty higher levels
        if len(levels) >= 6 and level >= 5:  # At least genus level available
            genus_level = levels[5] if len(levels) > 5 else ''
            species_level = levels[6] if len(levels) > 6 else ''
            
            # Extract genus name (remove g__ prefix if present)
            genus_name = genus_level[3:] if genus_level.startswith('g__') else genus_level
            species_name = species_level[3:] if species_level.startswith('s__') else species_level
            
            if level == 5:  # Genus level
                # For genus comparison, just return the genus name
                return genus_name if genus_name else unassignable_label
            elif level == 6:  # Species level
                # For species comparison, return genus + species
                if genus_name and species_name:
                    return f"{genus_name} {species_name}"
                elif genus_name:
                    return genus_name
                else:
                    return unassignable_label

        return ';'.join(levels[:level + 1])
    return f

def compute_taxon_accuracy(actual_table, expected_table, actual_sample_id=None,
                           expected_sample_id=None, ws_strip=True):
    """ Compute taxon accuracy and detection rates based on presence/absence of
    observations

        actual_table: table containing results achieved for query
        expected_table: table containing expected results
        actual_sample_id: sample_id to test (default is first sample id in
         actual_table.SampleIds)
        expected_sample_id: sample_id to test (default is first sample id in
         expected_table.SampleIds)
    """
    actual_obs_ids = get_observed_observation_ids(actual_table,
                                                  actual_sample_id,
                                                  ws_strip=ws_strip)
    expected_obs_ids = get_observed_observation_ids(expected_table,
                                                    expected_sample_id,
                                                    ws_strip=ws_strip)
    # length of list elements
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
    """ Return the set of observation ids with count > 0 in sample_id

        table: the biom table object to analyze
        sample_id: the sample_id to test (default is first sample id in
        table.SampleIds)
    """
    if sample_id is None:
        sample_id = table.ids(axis="sample")[0]

    result = []
    for observation_id in table.ids(axis="observation"):
        if table.get_value_by_ids(observation_id, sample_id) > 0.0:
            # remove all whitespace from observation_id
            if ws_strip is True:
                observation_id = "".join(observation_id.split())
            result.append(observation_id)

    return set(result)




def per_sequence_precision(expected_table_fp, actual_table_fp, feature_table,
                           sample_id, taxonomy_level, sample_weight=None, exclude=None):
    '''Precision/recall on individual representative sequences in a mock
    community.
    '''
    # locate expected and observed taxonomies
    exp_fp = expected_table_fp
    if exists(exp_fp):
        if exists(actual_table_fp):
            obs_fp = actual_table_fp
        else:
            raise RuntimeError('taxonomy assignments do not exist '
                               'for dataset {0}'.format(actual_table_fp))
        # compile lists of taxa only if observed in current sample
        exp = observations_to_list(exp_fp, feature_table, sample_id)
        obs = observations_to_list(obs_fp, feature_table, sample_id)
        try:
            exp, obs, exp_ids = _multiple_match_kludge(exp, obs, 
                                              fill_empty_observations=True)
        except AssertionError:
            print('AssertionError in:')
            print(actual_table_fp)
            raise
        # truncate taxonomies to the desired level
        exp_taxa, obs_taxa = load_prf(
            obs, exp, level=slice(0, taxonomy_level+1), sort=False)
        # compile sample weights (observations per sequence in sample)
        if sample_weight:
            sample_weight = [feature_table.get_value_by_ids(
                    line.split('\t')[0], sample_id) for line in exp]

        # run precision/recall
        total_exp, total_obs, ps, rs, fs, detailed_output = compute_prf(
            exp_taxa, obs_taxa, test_type='mock',
            sample_weight=sample_weight, exclude=exclude)
    else:
        total_exp, total_obs, ps, rs, fs, exp_ids, detailed_output = \
            -1., -1., -1., -1., -1., "None", "None"

    return total_exp, total_obs, ps, rs, fs, exp_ids, detailed_output


## Extract and import liens as a list
def observations_to_list(obs_fp, actual_table, sample_id):
    '''extract lines from obs_fp to list if they are observed in a given
    sample in biom table actual_table. Returns list of lines from obs_fp,
    which maps biom observation ids (first value in each line) to taxonomy
    labels in tab-delimited file.
    '''
    obs = import_to_list(input_data=obs_fp, field=slice(None, 2), delim="\t")
    obs = [line for line in obs if
           actual_table.exists(line.split('\t')[0], "observation") and
           actual_table.get_value_by_ids(line.split('\t')[0], sample_id) != 0]
    return obs

# Import lines of text file into list object
def import_to_list(input_data, field=None, delim='\t'):
    """Load lines from a file path or list, optionally extracting fields."""
    # Read from file if input is a string
    if isinstance(input_data, str):
        with open(input_data, "r") as f:
            lines = [line.strip() for line in f]
    elif isinstance(input_data, list):
        lines = [line.strip() for line in input_data]
    else:
        raise ValueError("input_data must be a file path or a list of lines")

    # Apply field extraction
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
    '''Sort expected and observed lists and kludge to deal with cases where we
    were unable to unambiguously select an expected taxonomy'''
    obs = {i: t for i, t in [r.split('\t', 1) for r in obs]}
    exp_grouped = defaultdict(list)
    for exp_id, exp_taxon in [r.split('\t') for r in exp]:
        exp_grouped[exp_id].append(exp_taxon)
    if fill_empty_observations:
        for k in exp_grouped.keys():
            if k not in obs.keys():
                obs[k] = 'Unassigned'
    else:
        assert obs.keys() == exp_grouped.keys(),\
            'observed and expected read labels differ:\n' + \
            str(list(obs.keys())) + '\n' + str(list(exp_grouped.keys()))
    new_exp = []
    new_obs = []
    exp_ids = []
    for exp_id, exp_taxons in exp_grouped.items():
        obs_row = '\t'.join([exp_id, obs[exp_id]])
        for exp_taxon in exp_taxons:
            # if obs[exp_id].startswith(exp_id):
            if obs[exp_id] == exp_taxon:
                row = '\t'.join([exp_id, exp_taxon])
                if len(exp_taxons) > 1:
                    print(f"\n Ambiguous expected taxonomy for ID: {exp_id}")
                    print(f" Expected:  {row}")
                    print(f" Observed:  {obs_row}")
                    print(" Candidate Expected taxonomies:")
                    for i, e in enumerate(exp_taxons, start=1):
                        print(f"  {i}. {e}")
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

    # raise error if obs_taxa and exp_taxa are not same length
    if len(obs_taxa) != len(exp_taxa):
        raise RuntimeError(
            'Lengths of expected and observed taxa do not match. '
            'Check inputs: {0}, {1}'.format(obs_fp, exp_fp))

    return exp_taxa, obs_taxa

def load_taxa(obs_fp, level=slice(0, 7), field=1, sort=False):
    '''Mount observed/expected taxonomy observations.
    obs_fp: path
        Input file containing taxonomy strings and IDs.
    level: slice
        Taxonomic range of interest. 0 = kingdom, 1 = phylum, 2 = class,
        3 = order, 4 = family, 5 = genus, 6 = species.
        Must use slice notation. For species, use level=slice(6, 7)
    field: int
        ab-delimited field containing taxonomy strings.
    '''
    raw_obs = import_to_list(obs_fp)
    if sort:
        raw_obs = sorted(raw_obs)
    obs = extract_taxa_names(raw_obs, field=field, level=level)
    obs = [';'.join([l.strip() for l in line.split(';')]) for line in obs
           if not line.startswith(('taxonomy', 'Taxon'))]
    return obs

def extract_taxa_names(infile, level=slice(6, 7), l_delim=';', field=None,
                       f_delim='\t', stripchars=None):
    '''Extract taxon names at a given level from taxonomy file OR LIST
    field = taxonomic level of interest. 0 = kingdom, 1 = phylum, 2 = class,
    3 = order, 4 = family, 5 = genus, 6 = species.
    Must use slice notation. For species, use level=slice(6, 7)

    stripchars: str
        Default, None, will strip leading and trailing whitespace. Set to ""
        to turn off character stripping.
    '''
    line_list = import_to_list(infile, field=field, delim=f_delim)

    # Truncate taxonomies and pass to set
    name_list = [l_delim.join(line.split(l_delim)[level]).strip(stripchars)
                 for line in line_list]
    return name_list

def compute_prf(exp, obs, test_type='cross-validated',
                l_range=range(1, 7), sample_weight=None, exclude=None):
    '''Compute precision, recall, and F-measure using sklearn.
    exp_taxa: list
        Expected observations for each sample (sequence).
    obs_taxa: list
        Predicted observations for each sample (sequence).
    avg: str
        Score averaging method using in sklearn. 'micro', 'weighted', or
        'macro'.
    test_type: str
        'novel-taxa' or 'cross-validated'.
    l_range: range
        Range of taxonomic levels to test if test_type = 'cross-validated'.
    sample_weight: array-like of shape = [n_samples], optional
        Sample weights.
    '''

    if test_type in ('mock', 'novel-taxa'):
        total_exp, total_obs, p, r, f, detailed = precision_recall_fscore(
            exp, obs, sample_weight=sample_weight, exclude=None)
    elif test_type == 'cross-validated':
        # initialize p/r/f as lists of 0s, representing each taxonomic level.
        p, r, f = [0] * 7, [0] * 7, [0] * 7
        # iterate over multiple taxonomic levels
        for lvl in l_range:
            _obs = extract_taxa_names(obs, level=slice(0, lvl+1))
            _exp = extract_taxa_names(exp, level=slice(0, lvl+1))
            p[lvl], r[lvl], f[lvl] = precision_recall_fscore(
                _exp, _obs, sample_weight=sample_weight, exclude=exclude)
    else:
        raise ValueError('test_type must be "novel-taxa" or "cross-validated" '
                         'or "mock".')

    return total_exp, total_obs, p, r, f, detailed

def precision_recall_fscore(exp, obs, sample_weight=None, exclude=None):
    """Compute precision, recall and F-measure with optional sample weights."""
    if exclude is None:
        exclude = []
    if sample_weight is None or sample_weight is False:
        sample_weight = [1] * len(exp)
    # length of list elements
    total_exp = len(exp)
    total_obs = len(obs)
    # detailed list
    detailed = []
    # compute score    
    tp, fp, fn = 0, 0, 0
    for e, o, w in zip(exp, obs, sample_weight):
        if e in exclude:
            continue
        classification = evaluate_classification(o, e)
        detailed.append(classification)
        if classification == 'match':
            # tp for the true class, the rest are tn
            tp += w
        elif classification == 'underclassification':
            # fn for the true class
            fn += w
            # no fp for the predicted class, because it was right to some level
            # the rest are tn
        else:
            # fp the the predicted class
            fp += w
            # fn for the true class, the rest are tn
            fn += w
    if tp == 0:
        return total_exp, total_obs, 0, 0, 0, ",".join(detailed)
    p = round((tp / (tp + fp)), 3)
    r = round((tp / (tp + fn)), 3)
    f = round((2*p*r / (p + r)), 3)
    return total_exp, total_obs, p, r, f, ",".join(detailed)

def evaluate_classification(obs_taxon, exp_taxon):
    '''Given an observed and actual taxonomy string corresponding to a cross-
    validated simulated community, score as match, overclassification,
    underclassification, or misclassification'''

    # Normalize IDTAXA's "unclassified_" format to empty levels for fair comparison
    def normalize_unclassified(taxon_string):
        """Convert IDTAXA's unclassified_prefix format to empty levels"""
        if not taxon_string:
            return taxon_string
            
        parts = [part.strip() for part in taxon_string.split(';')]
        normalized_parts = []
        
        for part in parts:
            if part.startswith('d__') or part.startswith('p__') or part.startswith('c__') or \
               part.startswith('o__') or part.startswith('f__') or part.startswith('g__') or \
               part.startswith('s__'):
                # Extract the content after the prefix
                prefix = part[:3]
                content = part[3:]
                
                # If content starts with "unclassified_", treat as empty level
                if content.startswith('unclassified_'):
                    normalized_parts.append(prefix)
                else:
                    normalized_parts.append(part)
            else:
                normalized_parts.append(part)
        
        return '; '.join(normalized_parts)
    
    # Normalize both taxonomies
    normalized_obs = normalize_unclassified(obs_taxon)
    normalized_exp = normalize_unclassified(exp_taxon)

    if normalized_obs == normalized_exp:
        return 'match'
    if normalized_exp.startswith(normalized_obs) or \
            normalized_obs in ('Unclassified', 'Unassigned', 'No blast hit'):
        return 'underclassification'
    if normalized_obs.startswith(normalized_exp):
        return 'overclassification'
    
    # Check for IDTAXA-style underclassification where observed has empty levels due to unclassified_
    obs_parts = [part.strip() for part in normalized_obs.split(';')]
    exp_parts = [part.strip() for part in normalized_exp.split(';')]
    
    if len(obs_parts) == len(exp_parts):
        # Check if observed taxonomy becomes "less specific" than expected after normalization
        # This happens when IDTAXA uses unclassified_ which gets normalized to empty levels
        obs_specificity = sum(1 for part in obs_parts if part and not part.endswith('__'))
        exp_specificity = sum(1 for part in exp_parts if part and not part.endswith('__'))
        
        if obs_specificity < exp_specificity:
            # Check if the non-empty parts match
            matches = 0
            for obs_part, exp_part in zip(obs_parts, exp_parts):
                if obs_part and not obs_part.endswith('__') and exp_part and not exp_part.endswith('__'):
                    if obs_part == exp_part:
                        matches += 1
                    else:
                        break  # Mismatch found
                elif obs_part and obs_part.endswith('__') and exp_part and not exp_part.endswith('__'):
                    break  # This is where classification stops in observed
            
            # If all matched parts are correct and observed is less specific
            if matches == obs_specificity:
                return 'underclassification'
    
    # Generic handling of genus/species classifications, treating ambiguous labels
    # such as "AMBIGUOUS" or empty strings as underclassification regardless of
    # the classifier that produced them
    obs_parts = [part.strip() for part in normalized_obs.split(';')]
    exp_parts = [part.strip() for part in normalized_exp.split(';')]

    if len(obs_parts) >= 6 and len(exp_parts) >= 6:
        ambiguous = {'', 'AMBIGUOUS', 'Unclassified', 'Unassigned', 'No blast hit'}

        obs_genus = obs_parts[5][3:] if obs_parts[5].startswith('g__') else obs_parts[5]
        exp_genus = exp_parts[5][3:] if exp_parts[5].startswith('g__') else exp_parts[5]

        # Missing or ambiguous genus is underclassification
        if obs_genus in ambiguous:
            return 'underclassification'

        if obs_genus == exp_genus:
            obs_has_species = len(obs_parts) > 6
            exp_has_species = len(exp_parts) > 6

            obs_species = obs_parts[6][3:] if obs_has_species and obs_parts[6].startswith('s__') else (obs_parts[6] if obs_has_species else '')
            exp_species = exp_parts[6][3:] if exp_has_species and exp_parts[6].startswith('s__') else (exp_parts[6] if exp_has_species else '')

            if exp_has_species:
                if obs_has_species:
                    if not obs_species or obs_species in ambiguous:
                        return 'underclassification'
                    return 'match' if obs_species == exp_species else 'misclassification'
                else:
                    return 'underclassification'
            else:  # expected genus level only
                if obs_has_species and obs_species and obs_species not in ambiguous:
                    return 'overclassification'
                return 'match'

    return 'misclassification'





def add_metadata_to_biom_table(input_fp, taxonomy_map_fp, biom_output_fp, mode='observed'):
    '''Load biom, add metadata, write to new table'''
    if mode == 'observed':
        newbiom = load_table(input_fp)
        if stat(taxonomy_map_fp).st_size == 0:
            metadata = {}
        else:
            metadata = MetadataMap.from_file(taxonomy_map_fp, header=['Feature ID', 'taxonomy', 'c'])
        newbiom.add_metadata(metadata, 'observation')
        write_biom_table(newbiom, 'json', biom_output_fp)

    elif mode == 'expected':
        # Load the abundance table amd process Taxon
        input_fp = pd.read_csv(input_fp, sep="\t")
        input_fp.set_index("Taxon", inplace=True)

        # prepare biom file
        newbiom = Table(
            input_fp.values,
            observation_ids=input_fp.index.tolist(),
            sample_ids=input_fp.columns.tolist(),
            observation_metadata=[{"taxonomy": idx} for idx in input_fp.index]
        )
        
        # write biom
        write_biom_table(newbiom, 'json', biom_output_fp)
        
        



def plot_from_data_frame(df, outdir, x_axis, y_vars, 
                         color_by, color_palette, line_by=None,
                         out_file="eval_results.pdf", style_theme="whitegrid"):
    '''Generate seaborn pointplot from pandas dataframe.
    df = pandas.DataFrame
    x_axis = x axis variable
    y_vars = LIST of variables to use for plotting y axis
    color_by = df variable on which to plot and color subgroups within data
    color_palette = color palette to use for plotting. Either a dict mapping
                     color_by groups to colors, or a named seaborn palette.
    style_theme = seaborn plot style theme
    plot_type = allows switching to other plot types, but this is untested
    '''
    # Create 2x3 subplot layout (adjust based on len(y_vars))
    ncols = 2
    nrows = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 12))
    axes = axes.flatten()

    sns.set_style(style_theme)
    sns.set_context("notebook", font_scale=1.2)

    # store handles/labels once
    legend_handles = None
    legend_labels = None


    for i, y_var in enumerate(y_vars):
        ax = axes[i]
        if line_by is None:
            sns.pointplot(
                data=df, x=x_axis, y=y_var, 
                hue=color_by, palette=color_palette, 
                scale=0.75, ax=ax
            )
        else:
            sns.lineplot(
                data=df, x=x_axis, y=y_var,
                hue=color_by, palette=color_palette,
                style=line_by,
                dashes=True,
                ax=ax
            )

        ax.set_title(y_var, fontsize=12)
        ax.set_xlabel(x_axis, fontsize=10)
        ax.set_ylabel(y_var, fontsize=10)
        ax.tick_params(labelsize=9)

        # Set consistent y-axis ticks and limits
        ax.set_ylim(0, 1.1)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

        # Save legend handles from the first plot only
        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

        # Remove individual subplot legends
        ax.legend_.remove()

    # Add shared legend in the bottom-right (unused) subplot
    if len(y_vars) < len(axes):
        legend_ax = axes[len(y_vars)]  # First unused axis
        legend_ax.axis('off')
        legend_ax.legend(
            legend_handles, legend_labels, title=color_by,
            fontsize=9, title_fontsize=10, loc='center', frameon=False
        )

    # Hide remaining unused subplots (if any)
    for j in range(len(y_vars) + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    # Save to PDF
    with PdfPages(join(outdir, out_file)) as pdf:
        pdf.savefig(fig)






def pointplot_from_data_frame_master(df, outdir, x_axis, y_vars, color_by,
                              color_palette, out_file="eval_results.pdf",
                              style_theme="whitegrid", plot_type=sns.pointplot):
    '''Generate seaborn pointplot from pandas dataframe.
    df = pandas.DataFrame
    x_axis = x axis variable
    y_vars = LIST of variables to use for plotting y axis
    color_by = df variable on which to plot and color subgroups within data
    color_palette = color palette to use for plotting. Either a dict mapping
                     color_by groups to colors, or a named seaborn palette.
    style_theme = seaborn plot style theme
    plot_type = allows switching to other plot types, but this is untested
    '''
    # Create 2x3 subplot layout (adjust based on len(y_vars))
    ncols = 2
    nrows = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, 12))
    axes = axes.flatten()

    sns.set_style(style_theme)
    sns.set_context("notebook", font_scale=1.2)

    # store handles/labels once
    legend_handles = None
    legend_labels = None

    for i, y_var in enumerate(y_vars):
        ax = axes[i]
        sns.pointplot(
            data=df, x=x_axis, y=y_var, hue=color_by,
            palette=color_palette, scale=0.75, ax=ax
        )
        ax.set_title(y_var, fontsize=12)
        ax.set_xlabel(x_axis, fontsize=10)
        ax.set_ylabel(y_var, fontsize=10)
        ax.tick_params(labelsize=9)

        # Set consistent y-axis ticks and limits
        ax.set_ylim(0, 1.1)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

        # Save legend handles from the first plot only
        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

        # Remove individual subplot legends
        ax.legend_.remove()

    # Add shared legend in the bottom-right (unused) subplot
    if len(y_vars) < len(axes):
        legend_ax = axes[len(y_vars)]  # First unused axis
        legend_ax.axis('off')
        legend_ax.legend(
            legend_handles, legend_labels, title=color_by,
            fontsize=9, title_fontsize=10, loc='center', frameon=False
        )

    # Hide remaining unused subplots (if any)
    for j in range(len(y_vars) + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()

    # Save to PDF
    with PdfPages(join(outdir, out_file)) as pdf:
        pdf.savefig(fig)