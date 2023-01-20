"""
module with generally useful functions used in multiple modules
"""

import sys

try:
    import requests
except BaseException:
    msg = "cannot import requests package, download_sample_abuns function not available"
    print(msg, file=sys.stderr)

from shutil import which


def mcl_available():
    """
    check whether MCL tool is in PATH and marked as an executable
    """
    return which("mcl") is not None


def download_sample_abuns(
        sample_id, out_fn,
        output_ids=['prot_ids'],
        prot_ids=[],
        id_type='prot_ids',
        url='https://www3.cmbi.umcn.nl/cedar/api/abundances'):
    """
    fetch abundances of a sample from CEDAR

    Args:
        sample_id (int): 
            CEDAR CRS number of a sample
        out_fn (str): 
            filepath, output location
        output_ids (list of strings, optional): Defaults to ['prot_ids'].
            the types of protein ids to include in the output.
            options: "prot_ids","prot_names","gene_names"
        prot_ids (list, optional): _description_. Defaults to [].
            if empty: fetches complete complexome profile
            otherwise: only fetch abundances for proteins matching given ids
        id_type (str, optional): Defaults to 'prot_ids'.
            type of identifier used when providing prot_ids
        url (str, optional): Defaults to 'https://www3.cmbi.umcn.nl/cedar/api/abundances'.
            url of CEDAR fetch_abundances api endpoint
    """
    headers = {
        'Content-Type': 'application/json',
        'Accept': 'text/csv',
    }
    data = {
        'sample_id': sample_id,
        'prot_ids': prot_ids,
        'id_type': id_type,
        'output_ids': output_ids,
    }
    response = requests.post(url, json=data, headers=headers)
    if response.ok:
        print(f'response ok, writing result to: {out_fn} ')
        with open(out_fn, 'wb') as f_obj:
            f_obj.write(response.content)


def map_df_index(df, mapping):
    """
    rename df's index using given mapping {index:new_id}

    original id is used for ids that have no mapping

    Args:
        df (pd.Dataframe): 
            table to map index of
        mapping (dict): 
            identifier mapping to use

    Returns:
        pd.Dataframe: table with mapped index
    """
    df = df.copy()
    df['mapped'] = df.index.map(mapping)

    # take care of missing labels
    missing = df['mapped'].isna()
    df.loc[missing, 'mapped'] = df.loc[missing].index.values

    df.set_index('mapped', inplace=True)
    return df


def get_stripped_mapping(full_id_list, sep='::'):
    """
    strip appendix from full ids

    stripped ids are stored in dict mapping back to full ids

    Args:
        full_id_list (list): 
            list of full ids to be stripped
        sep (str, optional): Defaults to '::'.
            separator between raw id and appendix

    Returns:
        dict: stripped ids as keys mapping to their original
              full ids with appendix
    """
    mapping = {}

    for full_id in full_id_list:
        stripped = full_id.rsplit(sep, 1)[0]
        if stripped in mapping.keys():
            mapping[stripped].append(full_id)
        else:
            mapping[stripped] = [full_id]
    return mapping


def invert_mapping(mapping_dict):
    """inverts given dict"""
    return {val: key for key, val in mapping_dict.items()}


def get_comp_mapping(left, right, nested_tags, mappings):
    """
    grab mapping for the current comparison of individual samples

    Args:
        left|right (str): 
            sample-level tags in this comparison
        nested_tags (dict): dict with nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        mappings (dict): contains id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject

    Returns:
        dict: mapping for given sample-level query,subject comparison
    """
    # get collection-level tags for current samples
    # OR: which collections do the samples belong to?
    # assumes it is present and only in 1 outer entry
    for key, values in nested_tags.items():
        if left in values:
            left_tag = key
        if right in values:
            right_tag = key

    # get (orthology) mapping if its provided
    mapping = get_col_mapping(left_tag, right_tag, mappings)

    return mapping


def get_col_mapping(left, right, mappings):
    """
    grab mapping for comparison of 2 collections from mappings

    Args:
        left/right (string): 
            collection identifiers
        mappings (dict of dicts): 
            dict with all available mappings

    Returns:
        dict: mapping between left and right collections
            keys: left collection ids
            values: corresponding right collection ids
    """
    if (left, right) in mappings.keys():
        mapping = mappings[(left, right)]
    elif (right, left) in mappings.keys():
        mapping = invert_mapping(mappings[(right, left)])
    else:
        mapping = False
    return mapping

def get_sample_tags(nested_tags):
    """
    get list of all replicate tags from nested_tags

    Args:
        nested_tags (dict of dicts): collection-replicate id structure

    Returns:
        list: all replicate-level ids
    """
    sample_tags = []
    for tags in nested_tags.values():
        sample_tags += tags
    return sample_tags

def get_comparison_matches(left, right, mapping=None):
    """
    determine id matches between comparison

    to get the id matches between indexes of
    to-be-compared samples, optionally using a mapping

    Args:
        left|right (list-like): 
            indexes of compared samples
        mapping (dict or None, optional): Defaults to None.
            dict with id mapping from left to right
            if None ids are directly compared

    Returns:
        list: ids that match between indexes
            OR
        dict: matching id pairs from left and right
            key: left id,  value: right id
    """
    if mapping:
        matches = {key: val for key, val in mapping.items()
                   if key in left and val in right}
    else:
        matches = [prot_id for prot_id in left if prot_id in right]

    return matches


def get_cluster_max_fraction(clusters, profile):
    """
    determine fraction in profile where cluster mean abundance is at max

    Args:
        clusters (dict): 
            cluster ids and lists of members
        profile (pd.DataFrame): 
            table with protein (rows) abundances
            in a number of fractions (columns)

    Returns:
        dict: for each cluster fraction at which mean abundance is at max value
    """
    # scale profile to weigh each protein equally
    scaled = profile.scale()

    maxfracs = {}
    for name, members in clusters.items():
        frac = scaled.loc[members].mean().reset_index(drop=True).idxmax()
        maxfracs[name] = frac

    return maxfracs


def correlate_samples(samples, method="pearson"):
    """
    compute interaction matrices for given samples

    Args:
        samples (dict of pd.df): samples to correlate
            samples are dataframes of feature data
            feature ids should be in index.
            column values should contain feature data
        method (str, optional): Defaults to "pearson".
            correlation method to be used.
            valid options: 'pearson','kendall','spearman'
    Returns:
        dict of pd.df: int_matrices
            contains symmetrical interaction matrix
            for each input sample
    """
    int_matrices = {}
    for name, sample in samples.items():
        correlated = sample.transpose().corr(method=method)
        int_matrices[name] = correlated

    return int_matrices


def eprint(*args, **kwargs):
    """
    like normal print but writes to stderror
    """
    print(*args, file=sys.stderr, **kwargs)
