import os
from os.path import join, basename, exists, dirname, splitext
from re import sub
import copy
from biom import load_table
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def delete_uncompleted_entries(database):
    """
    Delete entries without sequence or taxa information.

    Arguments
    ---------
        database: dict
            dictionary obtained from function `load_db_from_files()` or with the same format

    Returns
    -------
    dict
        Modified database dictionary

    See Also
    --------
    load_db_from_files : Load database into a dictionary
    download_db_from_dict : Download dictionary database

    """
    n_deleted = 0

    database_copy = copy.deepcopy(database)

    for ID, item in database.items():

        if len(item.keys()) != 2:
            # delete entry
            database_copy.pop(ID)

            # Info
            n_deleted = n_deleted + 1


    database = database_copy

    # Info
    if n_deleted > 0:
        print(f"{n_deleted} entries without sequence or taxa information were deleted")

    return database


def join_taxa_lineages(taxa):

    """
    Join list of clades into formatted taxa string, adding rank prefixes

    Arguments
    ---------
        taxa: list
            List of ranks to join

    Returns
    -------
    str

    See Also
    --------
    rm_prefs : Remove prefixes from a list of clades

    Examples
    -------
    >>> import utils
    >>> taxa = ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', 'Salmonella_enterica']
    >>> utils.join_taxa_lineages(taxa)
    'k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Salmonella; s__Salmonella_enterica'


    """

    prefixes = ["d__", "p__","c__","o__", "f__", "g__", "s__"]
    taxa_prefs = [''.join(item) for item in list(zip(prefixes, taxa))]
    taxa_full = '; '.join(taxa_prefs)
    return taxa_full


def rm_prefs(it, replace_unknowns = False):
    """
    Removes taxonomy prefixes from all the elements of the list

    Arguments
    ---------
        it: list or set or pandas dataframe column
            Taxonomy clades to modify
        replace_unknowns: bool, default = True
            If True, replaces undefined taxon levels with 'Unknown'.
            If False (default), those levels will become empty strings.

    Returns
    -------
    list
        Taxa without prefixes


    See Also
    --------
    join_taxa_lineages : join list of clades into formatted taxa string

    Examples
    --------
    >>> import utils
    >>> taxa = ['k__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacteriales', 'f__Enterobacteriaceae', 'g__Salmonella', 's__']
    >>> utils.rm_prefs(taxa)
    ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', '']

    >>> utils.rm_prefs(taxa, replace_unknowns = True)
    ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', 'Unknown']




    """
    result = [sub(r"^[dpcofgsd]__", "", item.strip().strip('"')) for item in it]

    if replace_unknowns:
        result = [i if i!= "" else "Unknown" for i in result]


    return result





def load_db_from_files(taxa_path=None, seqs_path=None):

    """
    Reads taxonomy file and sequence file from a database and stores it in
    a dictionary with the following format:

    .. code-block:: python

        {
            ID1: {'taxa': [rank1, rank2, rank3, ...], 'seq': Seq_object },
            ID2: ...
        }


    It also deletes entries with uncompleted informations (i.e. only taxa or seq)

    Arguments
    ---------
        taxa_path, seqs_path : str
            Paths of the taxa and seqs files of the database

    Returns
    -------
    dict
        Dictionary with database information


    See Also
    --------
    download_db_from_dict : Download dictionary database

    """

    def load_taxa(taxa_path, database):
        with open(taxa_path, "r") as taxafile:
            for line in taxafile:
                l_split = line.strip().split("\t")
                # taxa id
                ID = l_split[0].strip()

                if ID == 'Feature ID':
                    continue

                # Split Seven-lineages
                ftax = l_split[1].strip()
                ftax_split = rm_prefs(ftax.split(";"))

                database.setdefault(ID, {})
                database[ID]['taxa'] = ftax_split
        return database

    def load_seqs(seqs_path, database):
        with open(seqs_path, "r") as seqsfile:
            for record in SeqIO.parse(seqsfile, "fasta"):
                ID = record.id
                seq = record.seq

                database.setdefault(ID, {})
                database[ID]['seq'] = seq
        return database
 

    database = {}

    if taxa_path:
        database = load_taxa(taxa_path, database)
    if seqs_path:
        database = load_seqs(seqs_path, database)
    if taxa_path and seqs_path:
        database = delete_uncompleted_entries(database)

    return database



def download_db_from_dict(d, taxa_out_file=None, seqs_out_file=None):

    """
    Downloads database dictionary, creating a taxonomy file and/or a sequences file.
    The function will download only the provided files, e.g. if only the taxa file is provided,
    it will only download the database taxonomy.

    Arguments
    ---------
        d : dict
            database dictionary, as obtained by the function `load_db_from_files()`.
        taxa_out_file: str
            Taxonomy file output path
        seqs_out_file: str
            Sequence file output path

    Returns
    -------
    None

    Raises
    --------
    Exception
        If missing both file paths or missing the required taxa/sequence information in the database dictionary.


    See Also
    --------
    load_db_from_files : Load database files into a dictionary

    """


    records = []


    if taxa_out_file:
        with open(taxa_out_file, 'w') as taxa_out:
            # Add header to taxa file
            taxa_out.write("Feature ID\tTaxon\n")

            for ID, item in d.items():
                if taxa_out_file and 'taxa' in item:
                    taxa_full = join_taxa_lineages(taxa = item['taxa'])
                    taxa_out.write(f"{ID}\t{taxa_full}\n")
                elif taxa_out_file and 'taxa' not in item:
                    raise Exception(f"This dictionary doesn't contain taxa information. ID = {ID}")

                # Seqs file
                if seqs_out_file and 'seq' in item:
                    records.append(SeqRecord(item['seq'], id=ID, description=""))
                elif seqs_out_file and 'seq' not in item:
                    raise Exception(f"This dictionary doesn't contain sequence information. ID = {ID}")

    elif seqs_out_file:
        for ID, item in d.items():
            if seqs_out_file and 'seq' in item:
                    records.append(SeqRecord(item['seq'], id=ID, description=""))
            elif seqs_out_file and 'seq' not in item:
                raise Exception(f"This dictionary doesn't contain sequence information. ID = {ID}")


    else:
        raise Exception(f'No path provided.')

    # Write Seqs file
    if records:
        SeqIO.write(records, seqs_out_file,'fasta')



def check_dir(path):
    """
    Create ``path`` directories if they do not exist.

    Arguments
    ----------
        path : str
            Path of directories to check and create.

    Returns
    -------
    None

    """

    if not os.path.exists(path):
        os.makedirs(path)



def get_vregion_primers():
    """Get forward and reverse primers to extract 16S variable regions.

    Primers were extracted from this `article <https://doi.org/10.1128/msphere.01202-20>`_

    Return
    ------
    dict
        Dictionary with forward and reverse primers for 16S variable regions [1]_.

    References
    ----------
    .. [1] Abellan-Schneyder I, Matchado MS, Reitmeier S, Sommer A, Sewald Z, Baumbach J, List M, Neuhaus K. 2021. "Primer, pipelines, parameters: issues in 16S rRNA gene sequencing." mSphere 6:e01202-20., https://doi.org/10.1128/msphere.01202-20


    """
    vregion_dict = {
        'V1-V2': {"f":"AGAGTTTGATYMTGGCTCAG", "r": "GCTGCCTCCCGTAGGAGT"},
        'V1-V3': {"f":"AGAGTTTGATYMTGGCTCAG", "r": "ATTACCGCGGCTGCTGG"},
        'V3-V4': {"f":"CCTACGGGNGGCWGCAG", "r": "GACTACHVGGGTATCTAATCC"},
        'V3-V5': {"f": "CCTACGGGAGGCAGCAG", "r": "CCGTCAATTCMTTTRAGT"},
        'V4': {"f":"GTGCCAGCMGCCGCGGTAA", "r": "GGACTACNVGGGTWTCTAAT"},
        'V4-V5': {"f":"GTGCCAGCMGCCGCGGTAA", "r": "GAATTAAACCACATGCTC"},
        'V6-V8': {"f":"GAATTGACGGGGGCCCGCACAAG", "r": "CGGTGTGTACAAGGCCCGGGAACG"},
        'V7-V9': {"f":"CAACGAGCGCAACCCT", "r": "TACGGYTACCTTGTTACGACTT"},
        'V1-V9': {"f":"AGAGTTTGATCMTGGCTCAG", "r": "TACGGYTACCTTGTTACGACTT"},      
        'ITS1' : {"f": "CTTGGTCATTTAGAGGAAGTAA", "r": "GCTGCGTTCTTCATCGATGC"},
        'ITS2' : {"f": "GAACGCAGCRAANNGYGA", "r": "TCCTCCGCTTATTGATATGC"},
        'ITS12' : {"f": "CTTGGTCATTTAGAGGAAGTAA", "r": "TCCTCCGCTTATTGATATGC"},
    } 

    return vregion_dict



def get_vregion_lengths():
    """
    Get minimum and maximum length of the v regions

    Return
    ------
    dict
        Dictionary with the min and max length of V regions

    """
    vregion_dict = {
        'V1-V2': {'min': 310, 'max': 400},
        'V1-V3': {'min': 400, 'max': 600},
        'V3-V4': {'min': 400, 'max': 500},
        'V3-V5': {'min': 450, 'max': 600},
        'V4': {'min': 200, 'max': 350},
        'V4-V5': {'min': 400, 'max': 550},
        'V6-V8': {'min': 420, 'max': 600},
        'V7-V9': {'min': 450, 'max': 620},
        'V1-V9': {'min': 1300, 'max': 1550},
        'ITS1': {'min': 180, 'max': 350},
        'ITS2': {'min': 200, 'max': 400},
        'ITS12': {'min': 500, 'max': 800},
    }

    return vregion_dict
