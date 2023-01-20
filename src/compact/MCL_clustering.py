"""
create combined network from rbo scored data, cluster with MCL, and processes cluster results.

main functions
    - create_combined_network: given rbo scores between datasets, create combined network with normalised scores
    - run_MCL: runs MCL command line tool as subprocess
    - separate_subclusters: given MCL clusters, separates per collection, using replicates to compute fraction clustered scores
    - process_annot_MCL_res: processes MCL clusters into annotated CompaCt clusters and subclusters
"""

# base library imports
from itertools import combinations, chain, product
import subprocess as sp
from warnings import warn

# third party library
import numpy as np
import pandas as pd

# local library imports
from . import utils as ut
from . import process_data as prd
from .member_selection import select_members
from .utils import eprint,get_sample_tags

# check if mcl is available, warn if not
if not ut.mcl_available():
    warn("Warning. MCL not available, cannot perform"
         " clustering if MCL executable is not in PATH")

# normalise correlation and RBO scores

def normalise_scores(scores_list):
    """
    normalise top hit scores correcting for average per comparison

    Args:
        scores_list (list of pd series): top hits for each comparison

    Returns:
        list of pd series: normalised scores for each comparison
    """
    all_scores = pd.concat(scores_list)
    all_mean = all_scores.mean()

    comp_means = [comp.mean() for comp in scores_list]

    normalised = [scores / (comp_means[i] / all_mean)
                  for i, scores in enumerate(scores_list)]
    return normalised


def normalise_combined_scores(within_scores, between_scores, wbratio=1):
    """
    normalise within and between scores, average score per "comparison"

    Args:
        [within|between]_scores (list of pd series):
            top hit scores for each comparison
        wbratio (int, optional): Defaults to 1.
            ratio between average within and between scores
            determines ratio of average scores of within over between
            within-score-average will be wbratio*between-score-average

    Returns:
        tuple, containing normalised within and between scores
            scores: list of pd series
    """
    all_scores = pd.concat(within_scores + between_scores)
    all_mean = all_scores.mean()

    within_means = [comp.mean() for comp in within_scores]
    between_means = [comp.mean() for comp in between_scores]

    normalised_within = [wbratio * (scores / (within_means[i] / all_mean))
                         for i, scores in enumerate(within_scores)]
    normalised_between = [scores / (between_means[i] / all_mean)
                          for i, scores in enumerate(between_scores)]

    return normalised_within, normalised_between


def create_combined_network(between_scores, network_fn,
                            include_within=False,
                            within_scores=None,
                            wbratio=1):
    """
    generates network with combined edges from multiple comparisons

    Args:
        [within|between]_scores (list of pd series):
            top hit scores for each comparison.
        network_fn (string): 
            filepath of outut network
        include_within (bool, optional): Defaults to False.
            Whether to include within scores in combined network
        wbratio: numeric
            ratio between average within/between scores
            only used when include_within=True
    """
    # normalise the scores
    eprint('normalising comparison scores..')
    if include_within:
        norm_within, norm_between = normalise_combined_scores(
            within_scores, between_scores,
            wbratio=wbratio
        )

        # concatenate the comparisons' normalised scores
        total_network = pd.concat(norm_within + norm_between)
    else:
        norm_scores = normalise_scores(between_scores)
        total_network = pd.concat(norm_scores)

    # write total network to file (as abc file)
    total_network.to_csv(network_fn, sep='\t', header=False)

    eprint(f'combined network written to file: {network_fn}')

# perform MCL clustering using network with normalised edge weights


def run_MCL(input_fn, output_fn, inflation=2, processes=1):
    """
    runs mcl command line tool as subprocess

    Args:
        input_fn (string): 
            filepath of input network to be clustered
        output_fn (string):
            filepath of output result
        inflation (int, optional):
            mcl inflation param. Defaults to 2.
        processes (int, optional):
            number of processes/threads. Defaults to 1.
    """
    eprint('running MCL command line tool..')
    eprint(f'inflation parameter: {inflation}')
    eprint(f'input network file: {input_fn}')
    cmd = [
        'mcl',
        f'{input_fn}',
        '--abc',
        '-I',
        f'{inflation}',
        '-te',
        f'{processes}',
        '-o',
        f'{output_fn}']
    with sp.Popen(cmd, stdout=sp.PIPE, bufsize=1, universal_newlines=True) as p:
        for b in p.stdout:
            print(b, end='')
    eprint('MCL process ended\n')

# process the MCL results


def separate_subclusters(clusters, nested_tags):
    """
    separates clusters into subclusters per collection, pooling samples

    counts occurences of member in cluster for each sample

    Args:
        clusters (dict): all clusters
            keys: numeric cluster id
            value: list of cluster members
        nested_tags (dict of dicts): nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags

    Returns:
        dict of dicts: cluster members, separated per collection
            member occurence aggregated over samples per collection
    """
    clusters_split = {}
    for collection_id, tags in nested_tags.items():
        collection_clusters = {}
        for clust_id, clust_members in clusters.items():
            member_occurences = get_clust_member_occurence(clust_members, tags)
            if member_occurences:
                collection_clusters[clust_id] = member_occurences
        clusters_split[collection_id] = collection_clusters

    return clusters_split


def get_clust_member_occurence(clust, tags, as_fraction=True):
    """
    gets all collection members and occurence counts in given cluster

    Args:
        clust (list of strings): tagged cluster members
        tags (list of strings): sample tags to fetch and count proteins for
        as_fraction (bool, optional): Defaults to True.
            if True: returns occurence of members as fraction of total samples
            if False: returns occurence of members as raw counts

    Returns:
        dict: contains cluster members and their occurence
            key: members, value: occurence
    """
    all_members = []
    for tag in tags:
        repl_members = [prot[len(tag) + 1:] for prot in clust
                        if tag in prot]
        all_members += repl_members

    unique, counts = np.unique(all_members, return_counts=True)
    if as_fraction:
        fractions = counts / len(tags)
        return dict(zip(unique, fractions))
    return dict(zip(unique, counts))


def get_clust_info(clusters, clusters_split,
                   report_threshold=0.5):
    """
    creates table with cluster information

    Args:
        clusters (dict): all clusters
            keys: numeric cluster id
            value: list of cluster members
        clusters_split (dict of dicts): cluster members, separated per collection
            member occurence aggregated over samples per collection

        report_threshold (float, optional): Defaults to 0.5.
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold.

    Returns:
        pandas df: table with information per cluster
            rows: cluster ids
            columns: various cluster metrics
    """
    # initialize df
    clust_info = pd.DataFrame(index=clusters.keys())
    clust_info.index.name = 'clust_id'

    # add size of subclusters
    for tag, subclusters in clusters_split.items():
        sizes = {key: len(members) for key, members in
                 subclusters.items()}
        clust_info[f'{tag}_size'] = clust_info.index.map(sizes)

    # fill nas with zero
    clust_info.fillna(0, inplace=True)

    # add overall cluster sizes
    size_cols = [f'{tag}_size' for tag in clusters_split.keys()]
    clust_info['total_size'] = clust_info.loc[:, size_cols].sum(axis=1)

    # count number of collections represented
    presence = (clust_info.loc[:, size_cols] != 0)
    clust_info['n_represented'] = presence.sum(axis=1)

    # add subcluster sizes over report_threshold
    for tag, subclusters in clusters_split.items():
        threshold_sizes = {}
        for key, members in subclusters.items():
            threshold_members = [key for key, val in members.items()
                                 if val >= report_threshold]
            threshold_sizes[key] = len(threshold_members)
        colname = f'{tag}_over_{report_threshold}_size'
        clust_info[colname] = clust_info.index.map(threshold_sizes)

    # fill nas with zero
    clust_info.fillna(0, inplace=True)

    # add overall cluster sizes of proteins that make the threshold
    thresh_size_cols = [f'{tag}_over_{report_threshold}_size'
                        for tag in clusters_split.keys()]
    over_thresh_sizes = clust_info.loc[:, thresh_size_cols].sum(axis=1)
    clust_info[f'total_over_{report_threshold}_size'] = over_thresh_sizes

    # count number of collections represented over threshold
    thresh_presence = (clust_info.loc[:, thresh_size_cols] != 0)
    clust_info['robust_represented'] = thresh_presence.sum(axis=1)

    return clust_info


def invert_nested_tags(nested_tags):
    """
    creates sample to collection mapping from collection to sample

    Args:
        nested_tags (dict): collection-replicate mapping

    Returns:
        dict: replicate-collection mapping
    """
    inverted = {}
    for ctag, mtags in nested_tags.items():
        for mtag in mtags:
            inverted[mtag] = ctag
    return inverted


def fetch_tag_members(members, tag):
    """
    for list of tagged clust member ids,get stripped ids with given tag

    Args:
        members (list): cluser member ids, tagged with replicate ids
        tag (string): replicate id tag

    Returns:
        list: stripped ids of members that had given tag
    """
    tag_members = [mid[len(tag) + 1:] for mid in members
                   if tag in mid]
    return tag_members


def split_clusters(clusts, sample_tags):
    """
    split clusters into per-replicate subcluster members

    Args:
        clusts (dict): members for each cluster
        sample_tags (list): list with all replicate ids

    Returns:
        dict of dicts:
            keys: cluster id
            values: dict with member list for each replicate
    """
    nested_clusters = {}
    for cid, members in clusts.items():
        cur_clust = {}
        for tag in sample_tags:
            sample_members = fetch_tag_members(members, tag)
            if sample_members:
                cur_clust[tag] = sample_members
        nested_clusters[cid] = cur_clust
    return nested_clusters

def count_comp_matches(nested_cluster, mappings, nested_tags,
                       get_poss_matches=True):
    """
    count found matches per comparison for given cluster

    optionally also return the possible match counts per comparison

    Args:
        nested_cluster (dict): member list for each replicate
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        nested_tags (dict): collection-replicate id structure
        get_poss_matches (bool, optional): Defaults to True.
            whether to also return possible match counts per comparison

    Returns:
        dict: match counts for each comparison
    """
    sample_tags = get_sample_tags(nested_tags)
    comps = list(combinations(sample_tags, r=2))

    real_counts = {}
    poss_counts = {}
    for left, right in comps:
        # check if comp is present
        if not (left in nested_cluster.keys()
                and right in nested_cluster.keys()):
            continue

        left_members = prd.remove_appendices(nested_cluster[left])
        right_members = prd.remove_appendices(nested_cluster[right])
        comp_mapping = ut.get_comp_mapping(left, right, nested_tags, mappings)
        matches = ut.get_comparison_matches(
            left_members, right_members, mapping=comp_mapping
        )
        # get possible number of matches
        if get_poss_matches:
            n_poss = count_possible_matches(
                len(left_members), len(right_members))
            real_counts[(left, right)] = len(matches)
            poss_counts[(left, right)] = n_poss

    if get_poss_matches:
        return real_counts, poss_counts
    else:
        return real_counts


def count_possible_matches(n_left, n_right):
    """
    count possible matches for given comparison

    which is length of shortest list
    """
    if n_left < n_right:
        return n_left
    else:
        return n_right


def sum_counts(comp_counts, nested_tags):
    """
    sum per collection and total count

    in a way that loops over the total set of comparisons only once

    Args:
        comp_counts (dict): match counts for all comparisons
        nested_tags (dict): collection-replicate id structure

    Returns:
        dict: summed match counts for all collections as well as total
    """
    summed_counts = {col: 0 for col in nested_tags.keys()}
    summed_counts['total'] = 0
    sample_to_col = invert_nested_tags(nested_tags)
    for (left, right), count in comp_counts.items():
        left_col = sample_to_col[left]
        right_col = sample_to_col[right]

        # check if comparison is within collection
        if left_col == right_col:
            # if so add count to collection total
            summed_counts[left_col] += count

        # add count to total count
        summed_counts['total'] += count

    return summed_counts


def get_match_counts(nested_clusters, mappings, nested_tags,
                     report_fractions=True):
    """
    count matches between datasets/collections for every cluster

    Args:
        nested_clusters (dict of dicts):
            cluster members separated over replicates
                keys: cluster id
                values: dict with member list for each replicate

        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        nested_tags (dict): collection-replicate id structure
        report_fractions (bool, optional): Defaults to True.
            whether to only report raw counts or also normalised
            counts: divided by the number of possible matches
    Returns:
        pd.DataFrame: cluster match counts for all collections and total
            rows: clusters
            columns: collections and 'total'
    """
    match_counts = {}
    match_fractions = {}
    for cid, clust in nested_clusters.items():
        eprint(f'\rcounting matches for cluster: {cid}', end="")
        real_matches, poss_matches = count_comp_matches(
            clust, mappings, nested_tags)
        real_summed = sum_counts(real_matches, nested_tags)
        match_counts[cid] = real_summed
        # as fraction of possible matches
        if report_fractions:
            poss_summed = sum_counts(poss_matches, nested_tags)
            match_fractions[cid] = as_frac_of_poss(
                real_summed, poss_summed)

    match_counts = pd.DataFrame.from_dict(match_counts, orient='index')

    if report_fractions:
        match_fractions = pd.DataFrame.from_dict(
            match_fractions, orient='index')
        colnames = [f'{tag}_match_fraction' for tag in match_fractions.columns]
        match_fractions.columns = colnames
        return match_counts, match_fractions
    else:
        return match_counts


def as_frac_of_poss(real, poss):
    """
    convert counts to fraction of possible counts
    """
    frac_of_poss = {}
    for col, count in real.items():
        if count == 0:
            frac_of_poss[col] = 0
        else:
            frac_of_poss[col] = count / poss[col]
    return frac_of_poss


def fetch_subcluster_matches(subclusters, mappings, comps,
                             report_threshold=None):
    """
    find number of matches between subclusters

    Args:
        subclusters (dict): members for each subcluster of a cluster
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        comps (list): list containing comparison pairs
            elements (tuple): (left,right) sample tags to be compared
        report_threshold (float, optional): Defaults to 0.5.
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold.

    Returns:
        dict: list of matches for each comparison
            key (tuple): comparison with tags as (left,right)
            value (list of tuples): list of member matches: (left_id,right_id)
    """
    # check if comps are actually present in subclusters
    comps = [comp for comp in comps
             if comp[0] in subclusters.keys()
             and comp[1] in subclusters.keys()]

    comp_matches = {}
    # get mapping if provided
    for left, right in comps:

        # get mapping from stripped id to full ids for subclusters
        l_subclust = subclusters[left]
        r_subclust = subclusters[right]
        l_stripped_to_full = ut.get_stripped_mapping(l_subclust)
        r_stripped_to_full = ut.get_stripped_mapping(r_subclust)

        if (left, right) in mappings.keys():
            mapping = mappings[(left, right)]
        elif (right, left) in mappings.keys():
            mapping = prd.invert_mapping(mappings[(right, left)])
        else:
            mapping = False

        if mapping:
            # if l and r stripped ids are in mapping, add all
            # corresponding full id combinations to matches
            matches = []
            for l_stripped, l_full in l_stripped_to_full.items():
                if l_stripped in mapping.keys():
                    if mapping[l_stripped] in r_stripped_to_full.keys():
                        r_full = r_stripped_to_full[mapping[l_stripped]]
                        for l_id, r_id in product(l_full, r_full):
                            if report_threshold:  # in case of a report_threshold:
                                # skip if one of the proteins doesn't make the
                                # treshold
                                if (subclusters[left][l_id] < report_threshold
                                        or subclusters[right][r_id] < report_threshold):
                                    continue
                            pair = (f'{left}_{l_id}', f'{right}_{r_id}')
                            matches.append(pair)
        else:
            # if stripped l id is in stripped r ids, add all
            # corresponding full id combinations to matches
            matches = []
            for l_stripped, l_full in l_stripped_to_full.items():
                if l_stripped in r_stripped_to_full.keys():
                    r_full = r_stripped_to_full[l_stripped]
                    for l_id, r_id in product(l_full, r_full):
                        if report_threshold:
                            # skip if one of the proteins doesn't make the
                            # treshold
                            if (subclusters[left][l_id] < report_threshold
                                    or subclusters[right][r_id] < report_threshold):
                                continue
                        pair = (f'{left}_{l_id}', f'{right}_{r_id}')
                        matches.append(pair)

        comp_matches[(left, right)] = matches

    return comp_matches


def get_node_edge_tables(clusts, clusts_split, mappings, nested_tags, network):
    """
    generate node and edge tables for clustered nodes

    Args:
        clusts (dict): MCL output clusters
        clusts_split (nested dict):
                cluster members, split per collection, samples aggregated
                members are dict with key: member id, value: member weight
        network (pd.DataFrame): combined network stored as df of edges
            columns: left_node, right_node, edge weight

    Returns tuple with (nodes,edges):
        nodes: pd.Dataframe with annotated clustered nodes
        edges: pd.DataFrame with annotated edges between clustered nodes
    """
    # add cluster ids to edges, drop out of cluster edges
    eprint('adding clust ids to edges..')
    network = add_edge_clust_ids(network, clusts)

    # aggregate edges over samples
    eprint('aggregating edges over samples..')
    agg_network = aggregate_clust_edges(network, nested_tags)

    # add edge type column
    agg_network['edge_type'] = 'interaction'

    nodes = []

    comps = list(combinations(clusts_split.keys(), r=2))
    for i, clust_id in enumerate(clusts.keys()):
        eprint(
            f'processing cluster network {i+1} of {len(clusts)}..',
            end="\r")

        # get cluster nodes as dataframe, add to nodes
        node_df = get_cluster_nodes(clust_id, clusts_split)
        nodes.append(node_df)

        # get matches between subclusters
        subclusters = {}
        for tag in nested_tags.keys():
            if clust_id in clusts_split[tag].keys():
                subclusters[tag] = clusts_split[tag][clust_id]

        matches = fetch_subcluster_matches(subclusters, mappings, comps)

        # add matches as edges to aggregated network if there are any
        all_matches = list(chain.from_iterable(matches.values()))
        if any(True for _ in all_matches):
            all_matches = pd.DataFrame(all_matches)
            all_matches.columns = ['left', 'right']
            all_matches['edge_type'] = 'identity/ortholog'
            all_matches['clust_id'] = clust_id
            all_matches.set_index(['left', 'right'], inplace=True)
            agg_network = agg_network.append(all_matches)

    nodes = pd.concat(nodes)
    return nodes, agg_network


def select_clusters(clust_info, match_counts, min_match_count=2):
    """
    select clusters with members over threshold and at least two matches
    """
    robust_present = set(
        clust_info[clust_info['robust_represented'] != 0].index.values)
    nonzero_matches = set(
        match_counts[match_counts['total'] >= min_match_count].index.values)
    selected = list(robust_present & nonzero_matches)
    return selected


def process_annot_MCL_res(res_fn, nested_tags, network_fn,
                          mappings, report_threshold=0.5,
                          filter_clusters=True):
    """
    processes and annotates MCL results

    Parses raw MCL cluster results. Separates clusters per collection,
    aggregating membership over samples within a collection. Generates
    cluster info table, and optionally produces node and edge tables for
    network nodes and edges that are part of clusters for further analysis.

    Args:
        res_fn (str): 
            filepath of MCL output file
        nested_tags (dict of dicts): 
            nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        network_fn (str): 
            filepath, location of combined network file
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        report_threshold (float, optional): Defaults to 0.5.
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold.
        filter_clusters (bool, optional): default True
            whether to filter out clusters with less than 2 matches
            or proteins over report_threshold

    Returns:
        dict: processed mcl clustering results, containing:
            'clust_info': pd dataframe
                table with information of each identified cluster
            'clusts': dict
                clusters with all members together
            'clusts_split': (nested dict):
                cluster members, split per collection, samples aggregated
                members are dict with key: member id, value: member weight
            'best_guess': dict of dicts
                contains dict of subclusters for each collection
                subclusters: list of members that pass best guess selection
            'match_over_threshold': dict of dicts
                contains dict of subclusters for each collection
                subclusters: list of members that have match in another
                collection that is over fraction clustered threshold
            'edges': dataframe
                network edges part of one of the clusters
            'nodes': dataframe
                network nodes part of one of the clusters
    """
    eprint('processing MCL results..')
    # parse clusters
    clusts = prd.parse_MCL_result(res_fn)

    # split clusters per collection
    clusts_split = separate_subclusters(clusts,
                                        nested_tags)

    # get table with summary statistics of all clusters
    clust_info = get_clust_info(clusts, clusts_split,
                                report_threshold=report_threshold)

    # count matches between samples, total count and within collections
    sample_tags = []
    for tags in nested_tags.values():
        sample_tags += tags
    per_sample_clusters = split_clusters(clusts, sample_tags)
    match_counts, match_fractions = get_match_counts(
        per_sample_clusters, mappings, nested_tags)

    # filter clusters, must have matches and robustly present members
    if filter_clusters:
        selected = select_clusters(clust_info, match_counts)

        clusts = {key: val for key, val in clusts.items() if key in selected}
        new_clusts_split = {}
        for col, subclusts in clusts_split.items():
            new_subclusts = {key: val for key, val in subclusts.items()
                             if key in selected}
            new_clusts_split[col] = new_subclusts
        clusts_split = new_clusts_split
        clust_info = clust_info.loc[selected]
        match_fractions = match_fractions.loc[selected]

    # add match_fractions to clust_info
    clust_info = clust_info.merge(
        match_fractions,
        left_index=True,
        right_index=True)

    # get cluster node and edge tables for network analysis/visualization
    network = pd.read_csv(network_fn, sep='\t', header=None)
    nodes, edges = get_node_edge_tables(
        clusts, clusts_split, mappings, nested_tags, network)

    # get best guess members, and those that have an orth passing threshold
    best_guess, match_over_threshold = select_members(
        clusts_split, mappings,
        threshold=report_threshold
    )

    return {
        'clust_info': clust_info,
        'clusts': clusts,
        'clusts_split': clusts_split,
        'best_guess': best_guess,
        'match_over_threshold': match_over_threshold,
        'nodes': nodes,
        'edges': edges,
    }


def add_edge_clust_ids(network, clusts):
    """
    takes network table, adds clust_id annotation

    removes edges not associated with any cluster

    Args:
        network (pd df): containing network edges
            column structure: left_id,right_id,score
        clusts (dict): contains all MCL clusters

    Returns:
        pd df: network table with clust_id column added
    """
    network.columns = ['left', 'right', 'weight']
    network['clust_id'] = None

    for name, members in clusts.items():

        edge_in_clust = (network['left'].isin(members)
                         & network['right'].isin(members))
        network.loc[edge_in_clust, 'clust_id'] = name

    # drop edges not party of any cluster
    network.dropna(subset=['clust_id'], inplace=True)

    return network


def aggregate_clust_edges(edges, tags):
    """
    aggregates cluster edges over samples per collection

    Args:
        edges (pd df): network edge rows with clust ids
        tags (dict): nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags

    Returns:
        pd df: edges grouped by collection
            with average edge score and edge count
    """
    edges = edges.copy()

    # remove tags from ids
    # get ctags and rtags as separate cols
    for ctag, rep_tags in tags.items():
        for tag in rep_tags:
            to_strip = len(tag) + 1

            # left node col
            has_tag = edges['left'].str.startswith(tag)
            edges.loc[has_tag, 'left_rtag'] = tag
            edges.loc[has_tag, 'left'] = edges.loc[has_tag, 'left'].apply(
                lambda x: x[to_strip:])
            edges.loc[has_tag, 'left_ctag'] = ctag

            # right node col
            has_tag = edges['right'].str.startswith(tag)
            edges.loc[has_tag, 'right_rtag'] = tag
            edges.loc[has_tag, 'right'] = edges.loc[has_tag, 'right'].apply(
                lambda x: x[to_strip:])
            edges.loc[has_tag, 'right_ctag'] = ctag

    # aggregate edges over samples
    by = ['left', 'right', 'left_ctag', 'right_ctag', 'clust_id']
    grouped_edges = edges.groupby(by).weight.agg(
        ['mean', 'count']).reset_index()

    # add collection-level tags to node ids
    grouped_edges['left'] = grouped_edges['left_ctag'] + \
        '_' + grouped_edges['left']
    grouped_edges['right'] = grouped_edges['right_ctag'] + \
        '_' + grouped_edges['right']
    grouped_edges.drop(['left_ctag', 'right_ctag'], axis=1, inplace=True)

    # drop within sample edges
    rep_edges = (grouped_edges['left'] == grouped_edges['right'])
    grouped_edges = grouped_edges[~rep_edges]
    grouped_edges.set_index(['left', 'right'], inplace=True)

    return grouped_edges


def create_node_frame(nodes, ctag):
    """
    creates dataframe from cluster nodes

    Args:
        nodes (list of str): cluster members (nodes)
        ctag (str): collection level tag

    Returns:
        pd df: dataframe with nodes of current cluster
            columns: id, fraction_clustered, tag
    """
    as_df = pd.DataFrame(pd.Series(nodes)).reset_index()
    as_df.columns = ['id', 'fraction_clustered']
    tagged_id = as_df['id'].apply(lambda x: f'{ctag}_{x}')
    as_df.index = tagged_id
    as_df.index.name = 'tagged_id'
    as_df['tag'] = ctag
    return as_df


def get_cluster_nodes(clust_id, clusts_split):
    """
    for cluster with given clust_id, get node table

    Args:
        clust_id (int): id of current cluster
        clusts_split (dict): cluster members
            split per collection, samples aggregated

    Returns:
        pd df: node table
            columns: id,fraction_clustered,tag,clust_id
    """
    cluster_nodes = []
    for key, val in clusts_split.items():
        if clust_id in val.keys():
            cnodes = val[clust_id]
            as_frame = create_node_frame(cnodes, key)
            cluster_nodes.append(as_frame)
    node_frame = pd.concat(cluster_nodes)
    node_frame['clust_id'] = clust_id
    return node_frame
