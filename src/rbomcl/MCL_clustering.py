from itertools import combinations, chain, product

import pandas as pd
import numpy as np

import subprocess as sp

# local library imports
from . import process_data as prd
from .utils import get_stripped_mapping

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
                            include_within=True,
                            within_scores=None,
                            wbratio=1):
    """
    generates network with combined edges from multiple comparisons

    Args:
        [within|between]_scores (list of pd series):
            top hit scores for each comparison
        network_fn (string): filepath of outut network
        include_within (bool, optional): Defaults to True.
            Whether to include within scores in combined network
        wbratio: numeric
            ratio between average within/between scores
            only used when include_within=True
    """
    # normalise the scores
    print('normalising comparison scores..')
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

    print(f'combined network written to file: {network_fn}')

# perform MCL clustering using network with normalised edge weights

def run_MCL(input_fn, output_fn, inflation=2, processes=1):
    """
    runs mcl command line tool as subprocess

    Args:
        input_fn (string): filepath of input network to be clustered
        output_fn (string): filepath of output result
        inflation (int, optional): mcl inflation param. Defaults to 2.
        processes (int, optional): number of processes/threads. Defaults to 1.
    """
    print('running MCL command line tool..')
    print(f'inflation parameter: {inflation}')
    print(f'input network file: {input_fn}')
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
    print('MCL process ended\n')

# process the MCL results


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

    """
    returns df with cluster information

    (Kw)Args:
        pooled: whether result data has pooled samples
        report_threshold: at which threshold to count clust sizes
    """


def get_cluster_info(clusters, clusters_split,
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
    collection_presence = {}
    for i in clusters.keys():
        collection_presence[i] = {}
        for tag, cluster_part in clusters_split.items():
            if i in cluster_part.keys():
                collection_presence[i][tag] = True
            else:
                collection_presence[i][tag] = False

    cluster_info = pd.DataFrame.from_dict(
        collection_presence, orient='index')

    # add number of collections represented in cluster
    cluster_info['n_represented'] = cluster_info.sum(axis=1)

    # add size of collection parts
    for tag, subclusters in clusters_split.items():
        sizes = {key: len(members) for key, members in
                 subclusters.items()}
        cluster_info[f'{tag}_size'] = cluster_info.index.map(sizes)

        # add subcluster size above report_threshold
        threshold_sizes = {}
        for key, members in subclusters.items():
            threshold_members = [key for key, val in members.items()
                                 if val >= report_threshold]
            threshold_sizes[key] = len(threshold_members)
        colname = f'{tag}_over_{report_threshold}_size'
        cluster_info[colname] = cluster_info.index.map(threshold_sizes)

    # add overall cluster sizes
    size_cols = [f'{tag}_size' for tag in clusters_split.keys()]
    cluster_info['total_size'] = cluster_info.loc[:, size_cols].sum(axis=1)

    # add overall cluster sizes of proteins that make the threshold
    thresh_size_cols = [f'{tag}_over_{report_threshold}_size'
                        for tag in clusters_split.keys()]
    over_thresh_sizes = cluster_info.loc[:, thresh_size_cols].sum(axis=1)
    cluster_info[f'over_{report_threshold}_size'] = over_thresh_sizes

    return cluster_info


def process_MCL_result(res_fn, neted_tags,
                       report_threshold=0.5):
    """
    parses and processes raw MCL output

    Args:
        res_fn (string): filepath of MCL result
        nested_tags (dict of dicts): nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        report_threshold (float, optional): Defaults to 0.5.
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold.

    Returns:
        clusts (dict): contains all MCL clusters
        clusts_split (dict): cluster members
            split per collection, samples aggregated
        clust_info (pd df): table with information per cluster
            rows: cluster ids
            columns: various cluster metrics./=
    """
    clusts = parse_MCL_result(res_fn)
    clusts_split = separate_subclusters(clusts,
                                        neted_tags)
    # STILL NEED TO REWRITE separates TOOL TO PROPERLY COUNT
    # TOTAL CLUSTER SIZES (AFTER POOLING)
    clust_info = get_cluster_info(clusts, clusts_split,
                                  report_threshold=report_threshold)
    return clusts, clusts_split, clust_info

    """
    find number of matches between subclusters

    subclusters: dict with tag:subcluster structure
    """


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
        l_stripped_to_full = get_stripped_mapping(l_subclust)
        r_stripped_to_full = get_stripped_mapping(r_subclust)

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


def avg_clust_weight(size, total_weight):
    """
    computes average edge weight of cluster

    Args:
        size (int): total cluster size
        total_weight (float): total edge weight

    Returns:
        float: average cluster edge weight
            as total_edge_weight/n_possible_edges
    """
    if size == 1:
        return 0
    poss_edges = (size * size - size) / 2
    return total_weight / poss_edges


def get_match_counts(matches):
    """
    counts number of matches for each comparison

    Args:
        matches (dict): matches per comparison

    Returns:
        dict: number of matches for each comparison
    """
    return {comp: len(comp_matches) for comp, comp_matches in matches.items()}

    """
    version of annotate_MCL results to work for pooled results
    """


def process_annot_MCL_res(res_fn, nested_tags, network_fn,
                          mappings, report_threshold=0.5):
    """
    processes and annotates MCL results

    Parses raw MCL cluster results. Separates clusters per collection,
    aggregating membership over samples within a collection. Generates
    cluster info table, and produces node and edge tables for network
    nodes and edges that are part of clusters for further analysis.

    Args:
        res_fn (str): filepath of MCL output file
        nested_tags (dict of dicts): nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        network_fn (str): filepath, location of combined network file
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        report_threshold (float, optional): Defaults to 0.5.
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold.

    Returns:
        dict: processed mcl clustering results, containing:
            'clust_info': pd dataframe
                table with information of each identified cluster
            'clusts': dict
                clusters with all members together
            'clusts_split': dict
                cluster members, split per collection, samples aggregated
            'edges': dataframe
                network edges part of one of the clusters
            'nodes': dataframe
                network nodes part of one of the clusters
    """
    print('processing MCL results..')
    clusts, clusts_split, clust_info = process_MCL_result(
        res_fn, nested_tags, report_threshold=report_threshold)
    network = pd.read_csv(network_fn, sep='\t', header=None)

    # add cluster ids to edges, drop out of cluster edges
    print('adding clust ids to edges..')
    network = add_edge_clust_ids(network, clusts)

    # aggregate edges over samples
    print('aggregating edges over samples..')
    agg_network = aggregate_clust_edges(network, nested_tags)

    # add edge type column
    agg_network['edge_type'] = 'interaction'

    nodes = []

    comps = list(combinations(clusts_split.keys(), r=2))
    for clust_id in clusts.keys():
        print(f'processing cluster {clust_id+1} of {len(clusts)}..', end="\r")

        # compute average edge weight
        size = len(clusts[clust_id])
        total_weight = network.loc[network['clust_id']
                                   == clust_id, 'weight'].sum()
        avg_weight = avg_clust_weight(size, total_weight)
        clust_info.loc[clust_id, 'avg_weight'] = avg_weight

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

        thresh_matches = fetch_subcluster_matches(
            subclusters, mappings, comps, report_threshold=report_threshold)

        # add match counts per comparison to clust info
        match_counts = get_match_counts(matches)
        total_matches = sum(match_counts.values())
        for comp, count in match_counts.items():
            clust_info.loc[clust_id, f'{comp[0]}:{comp[1]}_matches'] = count
        clust_info.loc[clust_id, 'total_matches'] = total_matches

        # add over threshold matches to clust info
        thresh_match_counts = get_match_counts(thresh_matches)
        total_thresh_matches = sum(thresh_match_counts.values())
        for comp, count in thresh_match_counts.items():
            clust_info.loc[clust_id,
                           f'{comp[0]}:{comp[1]}_over_{report_threshold}_matches'] = count
        clust_info.loc[clust_id,
                       f'total_over_{report_threshold}_matches'] = total_thresh_matches

        # add threshold based stringent presence annotation
        filter_subcluster_presence(clust_info, nested_tags.keys(),
                                   report_threshold)

    return {
        'clust_info': clust_info,
        'clusts': clusts,
        'clusts_split': clusts_split,
        'edges': agg_network,
        'nodes': pd.concat(nodes)
    }


def filter_subcluster_presence(clust_info, tags, report_threshold):
    """
    add threshold-based subcluster annotations to clust_info table

    Args:
        clust_info (pd df): table with information of each identified cluster
        tags (list-like): collection-level tags
        report_threshold (float, optional): Defaults to 0.5.
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold.
    """
    for tag in tags:
        # are they represented over threshold?
        clust_info[f'{tag}_present'] = clust_info[f'{tag}_over_{report_threshold}_size'] > 0

        # are there matches with each collection over threshold?
        match_cols = [col for col in clust_info.columns
                      if tag in col and f'over_{report_threshold}_matches' in col]

        clust_info[f'{tag}_has_match'] = clust_info[match_cols].sum(axis=1) > 0

    # how many complecomes are "robustly" present
    has_match_cols = [f'{tag}_has_match' for tag in tags]
    clust_info['robust_n_represented'] = clust_info[has_match_cols].sum(axis=1)


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
            columns: id, fraction_present, tag
    """
    as_df = pd.DataFrame(pd.Series(nodes)).reset_index()
    as_df.columns = ['id', 'fraction_present']
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
            columns: id,fraction_present,tag,clust_id
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


if __name__ == "__main__":

    pooled_tags = {
        'GAM': ['GAM1', 'GAM2', 'GAM3', 'GAM4'],
        'AS': ['AS1', 'AS2', 'AS3', 'AS4'],
        'TG': ['TG']
    }

    # # annot_network = pd.read_csv('../data/TEST_network_edges.tsv',sep='\t')
    # print(annot_network.head())

    # aggregate_clust_edges(annot_network,pooled_tags)

    # quit()

    mclres_fn = '../results/GAM_AS_TG_separated_reps_results_laptop/mcl_result.tsv'

    network_fn = '../results/GAM_AS_TG_separated_reps_results_laptop/combined_network.tsv'

    orths = prd.parse_orthology('../data/PF3D7_to_TGGT1_consensus_copied.tsv')
    mappings = {
        ('AS', 'TG'): orths,
        ('GAM', 'TG'): orths,
    }

    network_annot = process_annot_MCL_res(
        mclres_fn, pooled_tags, network_fn, mappings)

    print(network_annot.keys())
    # print(network_annot.head())
    # network_annot.to_csv('../data/TEST_network_edges.tsv',sep='\t',index=False)

    quit()
    # edges = pd.read_csv('../data/TEST_clust_edges.tsv',sep='\t')
    # edges.columns = ['left','right','weight']

    # res = aggregate_clust_edges(edges,pooled_tags)
    # print(res)

    # quit()

    clust_data = res['clust_data']

    # aggregate_pooled_edges(clust_data,pooled_tags)

    cur_clust = 26
    # clust_edges = clust_data[cur_clust]['edges']
    # clust_edges.to_csv('../data/TEST_clust_edges.tsv',sep='\t', index=False)

    """
    # NOW FOR MERCEDES' clusters
    import process_eupath_orthology as po

    mcl_result_fn = '../results/merc_pooled_3_species_results/mcl_result.tsv'
    combined_net_fn = '../results/merc_pooled_3_species_results/combined_network.tsv'
    ANKA_samples = [f'CRS{i}' for i in range(58,64)]
    PKNH_samples = [f'CRS{i}' for i in range(64,70)]
    PF3D7_samples = [f'CRS{i}' for i in range(70,76)]

    nested_tags = {
        'PF3':PF3D7_samples,
        'BER':ANKA_samples,
        'KNO':PKNH_samples,
    }

    PF3D7_ANKA = po.process_orths('../data/PF3D7_to_ANKA_orthoMCL.tsv',species_tag='PBANKA')
    PF3D7_PKNH = po.process_orths('../data/PF3D7_to_PKNH_orthoMCL.tsv',species_tag='PKNH')
    PKNH_ANKA = po.process_orths('../data/PKNH_to_ANKA_orthoMCL.tsv',species_tag='PBANKA')

    mappings = {
        ('PF3','BER'):PF3D7_ANKA,
        ('PF3','KNO'):PF3D7_PKNH,
        ('KNO','BER'):PKNH_ANKA,
    }

    res = process_annot_MCL_res(mcl_result_fn,nested_tags,combined_net_fn,
                                     mappings)

    cur_clust = 17
    print(res['clust_data'][cur_clust]['edges'])
    print()
    print(res['clusts'][cur_clust])
    """