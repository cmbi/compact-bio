"""
general functions for processing,parsing data etc.
"""

from collections import Counter

import pandas as pd
import numpy as np

TOXO_ANNOT_FN = '../data/TGGT1_annotation.tsv'
PLASMO_ANNOT_FN = '../data/PF3D7_annotation.tsv'
PBANKA_ANNOT_FN = '../data/PBANKA_annotation.tsv'
PKNH_ANNOT_FN = '../data/PKNH_annotation.tsv'


def parse_profile(tsv_fn):
    """
    parses standard format complexome profile into df

    Args:
        tsv_fn (str): filepath of file containing:
            complexome profile in tab separated text format.
            single header col with fraction ids
            single index row with protein ids
            numeric abundance values
    Returns:
        pd df: complexome profile as dataframe
    """
    df = pd.read_csv(tsv_fn, sep='\t')
    df.iloc[:, 0] = df.iloc[:, 0].astype(str)
    df.set_index(df.columns[0], inplace=True)
    df.index.name = 'prot_ids'
    return df


def parse_mapping(tsv_fn):
    """
    parse file with identifier mapping into dict

    Args:
        tsv_fn (str): filepath of mapping file
            each line contains 2 tab-separated ids
    Returns:
        dict: id mapping as  {left:right,..}
    """
    with open(tsv_fn, 'r') as f_obj:
        as_dict = {}
        for line in f_obj:
            key, val = line.strip().split('\t')
            as_dict[key] = val
        return as_dict


def fetch_mapping(left_tag, right_tag, mappings):
    """
    fetches correct mapping for left:right from mappings

    inverts mapping if only right:left mapping is present
    returns None if no mapping is available

    Args:
        [left|right]_tag (str): collection level tags
            for which to fetch mapping
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject

    Returns:
        dict: dict with mapping for given ids
    """
    if (left_tag, right_tag) in mappings.keys():
        mapping = mappings[(left_tag, right_tag)]
    elif (right_tag, left_tag) in mappings.keys():
        mapping = invert_mapping(mappings[(right_tag, left_tag)])
    else:
        mapping = None
    return mapping


def invert_mapping(mapping_dict):
    """inverts given dict"""
    return {val: key for key, val in mapping_dict.items()}


def split_prot_ids(df):
    """
    if index contains rows with multiple protein ids, split them

    Args:
        df (pd df): protein df for which index should be split

    Returns:
        list of lists: list of row indices, each index is list
        with one or more ids
    """
    ids = df.reset_index()['prot_ids'].apply(lambda x: x.split(','))
    return ids.to_list()


def match_ids(ids_list, target_list):
    """
    match list with given ids to target list

    to be used when dealing with data that have multiple ids
    per row that need to be narrowed down to 1. It will pick
    one of the ids per row that ensures the most matches with
    ids in target list
    CARE: IF MULTIPLE MATCHES WITH TARGET, ARBITRARILY PICKS FIRST

    Args
        ids_list (list-of-lists):
        target_list (list): list of ids
    Returns
        matched_list: list

    """
    matched_list = []
    for ids in ids_list:

        if len(ids) == 0:
            msg = 'row without id?!'
            raise ValueError(msg)
        matches = [i for i in ids if i in target_list]
        if len(matches) == 1:
            matched_list.append(matches[0])
        elif len(matches) == 0:
            matched_list.append(ids[0])
        else:
            matched_list.append(matches[0])
    return matched_list


def split_match_ids(df, target_list):
    """
    reduces df with multiple ids per row to 1 id per row

    picks id so that matches with target_list is maximized

    Args:
        df (pd df): index that can have multiple
            ids per row, in comma-separated strings
        target_list (list): list of ids to which
            dataframe index should be matched

    Returns:
        pd df: copy of input df with single id per row
    """
    split_ids = split_prot_ids(df)
    matched_singles = match_ids(split_ids, target_list)
    copy = df.copy()
    copy.index = matched_singles
    return copy


def parse_orthology(orth_fn):
    df = pd.read_csv(orth_fn, sep='\t', index_col=0, header=None)
    return df[1].to_dict()


def parse_annotation(annot_fn):
    return pd.read_csv(annot_fn, sep='\t', index_col=0)


def parse_scores(score_fn, rename_dups=True):
    """
    parses matrix-structured scores into dataframe

    Args:
        score_fn: path/fn of tsv table with scores
                  matrix structure with index and columns

    Returns: pd.DataFrame
    """
    df = pd.read_csv(score_fn, sep='\t', index_col=0)

    # rename duplicate occurences to get unique indices
    if rename_dups:
        df.index = rename_duplicates(list(df.index.values))
        df.columns = rename_duplicates(list(df.columns))

    return df


def parse_top_hits(top_hit_fn):
    return pd.read_csv(
        top_hit_fn,
        sep='\t',
        header=None,
        index_col=[0, 1],
        squeeze=True
    )

# FUNCTIONALITY AROUND COMBINED NETWORK


def parse_network(net_fn):
    return pd.read_csv(net_fn, sep='\t', header=None)


def write_network(net, out_fn):
    net.to_csv(out_fn, sep='\t', index=False, header=False)


def filter_network(network, comps):
    """
    filter provided comparisons from given network

    Args:
        network (pd df): network containing multiple comparisons
        comps (list of tuples):
            comparisons to include in the output network
            comparison (tuple, (str,str)):
                left and right tags of compared samples

    Returns:
        pd df: subnetwork containing only provided comparisons
    """
    comp_data = []
    for left, right in comps:
        is_comp = (network[0].str.startswith(left)
                   & network[1].str.startswith(right))
        # assuming comparison is either stored completely
        # left --> right or completely right --> left
        if is_comp.sum() == 0:
            is_comp = (network[0].str.startswith(right)
                       & network[1].str.startswith(left))

        # if reverse also yields nothing, comparison not present!
        if is_comp.sum() == 0:
            msg = f'comparison ("{left}":"{right}") not present in input network'
            raise ValueError(msg)

        comp_data.append(network.loc[is_comp])

    sub_net = pd.concat(comp_data)
    return sub_net


def parse_MCL_result(res_fn):
    """
    parses MCL result from file into dict

    Args:
        res_fn (string): filepath of MCL result

    Returns:
        dict: containing all MCL clusters
    """
    clusters = {}
    with open(res_fn, 'r') as f_obj:
        for i, line in enumerate(f_obj):
            ids = line.strip().split('\t')
            # ignore clusters with only one element
            if len(ids) > 1:
                clusters[i] = ids
    return clusters


def annotate_df(to_annot, annot_fn):
    """
    """
    annotation = parse_annotation(annot_fn)
    merged = to_annot.merge(
        annotation,
        how='left',
        left_index=True,
        right_index=True)
    return merged


def rename_duplicates(sorted_list, separator="::"):
    """
    in sorted list of ids, renames occurences after first

    Args:
        sorted_list (list): list of string identifiers
        separator (str, optional): Defaults to "::".
            str that will separate original id and
            number that will be appended

    Returns:
        list: list with duplicate ids renamed
    """
    new_list = sorted_list.copy()
    # find not-unique protein ids
    not_unique = [prot for prot, occ in
                  Counter(sorted_list).items() if occ > 1]
    for prot in not_unique:
        # find locations of not-unique protein
        indices = [i for i, x in enumerate(sorted_list)
                   if x == prot]

        # add suffix to occurences after first
        for i, idx in enumerate(indices[1:]):
            suffix = f'{separator}{i+1}'
            new_list[idx] = prot + suffix

    return new_list


def remove_appendices(id_list, separator="::"):
    """
    removes trailing appendices from ids

    Args:
        id_list (list): ids that might have appendix
        separator (str, optional): Defaults to "::".
            everything including and after separator will
            be stripped from the string

    Returns:
        list of strings: stripped ids
    """
    return [name.rsplit(separator, 1)[0] for name in id_list]


def rename_duplicates_int_matrix(df, separator="::"):
    """
    renames duplicate row and col ids inplace in currenf df

    Args:
        df (pd df): dataframe from which duplicate
            column and rows will be renamed

        separator (str, optional): Defaults to "::".
            everything including and after separator will
            be stripped from the string

    Returns:
        int: number of duplicate ids in df index
    """
    n_dups = df.index.duplicated(keep=False).sum()

    # rename duplicate occurences to get unique indices
    df.index = rename_duplicates(list(df.index.values),
                                 separator=separator)
    df.columns = rename_duplicates(list(df.columns),
                                   separator=separator)

    return n_dups


def add_tag_df(rowtag, coltag, matrix):
    """
    prepends the collection tag to each protein id

    applied to both rows and columns of the matrix

    Args:
        [row|col]tag (str): tag to prepend to ids
        matrix (pd df): matrix with index and columns

    Returns:
        pd df: copy of input df with tagged ids
    """

    # create vectorized function to add tags to array
    add_rowtag = np.vectorize(lambda x: f'{rowtag}_{x}')
    add_coltag = np.vectorize(lambda x: f'{coltag}_{x}')

    new_index = add_rowtag(matrix.index.values)
    new_columns = add_coltag(matrix.columns.values)

    new_matrix = matrix.copy()

    new_matrix.index = new_index
    new_matrix.columns = new_columns

    return new_matrix


def add_tag_multiindex(left_tag, right_tag, multiindex):
    """
    prepend complexome tag to protein ids in 2-level multiindex

    Args:
        [left|right]_tag (str): tag to prepend to ids
        multiindex (pd MultiIndex): 2-level multiindex,
            tags will be prepended to both levels

    Returns:
        pd MultiIndex: new index with tagged ids
    """
    add_left_tag = np.vectorize(lambda x: f'{left_tag}_{x}')
    add_right_tag = np.vectorize(lambda x: f'{right_tag}_{x}')

    multiindex = multiindex.set_levels(
        add_left_tag(multiindex.levels[0]),
        level=0)
    multiindex = multiindex.set_levels(
        add_right_tag(multiindex.levels[1]),
        level=1)

    return multiindex


def split_clustmember_tables(nodes, mappings):
    """
    get complexome-specific cluster membership tables from nodes

    uses mapping to add mappings to other complexomes

    Args:
        nodes (pd df): table of cluster nodes
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject

    Returns:
        dict: containing separate table of nodes per tag
            key (str): tag, value (pd df): node table
    """
    # set clust_ids to index, split tables per tag
    nodes = nodes.set_index('clust_id')

    # split tables up into different dataframes
    split_tables = {}
    tags = nodes['tag'].unique()
    for tag in tags:
        frame = nodes.loc[nodes['tag'] == tag].drop('tag', axis=1)

        # add mappings to other complexomes
        for match in tags:
            if match == tag:
                continue
            # add mapping ids if pair has mapping
            mapping = fetch_mapping(tag, match, mappings)
            if mapping:
                frame[f'{match}_mapping'] = frame['id'].map(mapping)
        split_tables[tag] = frame

    return split_tables


def parse_gmt(filename):
    """
    parses gmt format file with named reference groups

    Args:
        filename (string): filepath of .gmt file

    Returns:
        dict: dict with group names and members
    """
    with open(filename, 'r', errors='ignore') as file_object:
        complex_dict = dict()
        for line in file_object:
            linelist = line.strip().split('\t')
            name = linelist[0]
            members = linelist[2:]
            complex_dict[name] = members

    return complex_dict


if __name__ == "__main__":

    # parse combined network
    network_fn = '/home/joerivs/Documents/Apicomplexa_project/results/c12_run_gene_Apr28_results/test_chunk_network.tsv'
    network_fn = '/home/joerivs/Documents/Apicomplexa_project/results/c12_run_gene_Apr28_results/combined_network.tsv'

    net = parse_network(network_fn)
    print(net.head())

    # filter combined network

    comps = [('CRS86', 'CRS50')]
    subnet = filter_network(net, comps)

    print(subnet.shape)
    print(subnet.head())

    # write combined network
    write_network(
        subnet,
        '/home/joerivs/Documents/Apicomplexa_project/results/c12_run_gene_Apr28_results/CRS86_CRS50_subnet.tsv')

    # extract subnetwork

    # parse dataframe to annotate, for testing
    # to_annot_fn = '../results/top_03_percent_top_hit_clusters_PLASMO.tsv'
    # to_annot = pd.read_csv(to_annot_fn, sep='\t', index_col=0)
    # to_annot = to_annot[['average_rbo']]
    # annotate_df(to_annot, PLASMO_ANNOT_FN)
