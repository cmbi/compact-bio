"""
This module contains the high-level workflow of this project,
to run the major steps of the CompaCt analsysis.

The main functions:
    - between_scoring: rbo scoring between all provided correlation datasets
    - score_comparison: rbo scoring and top hit selection between two correlation datasete
    - mcl_clustering: perform clustering on computed rbo scores
    - process_mcl_result: process raw MCL output to get annotated CompaCt result clusters and subclusters
    - save_results: save CompaCt results to an output folder of choice
"""

# base imports
from multiprocessing.sharedctypes import Value
import os
from itertools import combinations

# local library imports
from .cluster_annotation import annotate_clusters
from . import process_data as prd
from . import pairwise_scoring as ps
from . import reciprocal_hits as rh
from . import MCL_clustering as mcl
from . import utils as ut
from .utils import eprint

### Helper Functions ###


def report_matches(comps, int_matrices, nested_tags, mappings, min_count=10):
    """
    reports low match numbers for each sample comparison

    Args:
        comps (list): list containing comparison pairs
            elements (tuple): (left,right) sample tags to be compared
        int_matrices (dict): contains interaction matrices
            keys: sample-level tag of interaction matrix
            values: pd dataframe with pairwise interaction scores
        nested_tags (dict): nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        mappings (dict): id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        min_count (int, optional): Defaults to 10.
            match numbers lower than this value will be reported

    Raises:
        ValueError: in case of invalid min_count parameter
    """
    for left, right in comps:
        mapping = ut.get_comp_mapping(left, right, nested_tags, mappings)

        left_index = int_matrices[left].index
        right_index = int_matrices[right].index

        matches = ut.get_comparison_matches(left_index, right_index,
                                            mapping=mapping)
        match_count = len(matches)
        if min_count is None:
            eprint(
                f'number of matches between {left} and {right}: {match_count}')
        elif match_count < min_count:
            eprint(
                f'warning! low number of matches between {left} and {right}: {match_count}')
        elif match_count >= min_count:
            pass
        else:
            raise ValueError(
                f"min_count must be None or int, not: {min_count}")


def check_matches(nested_tags, int_matrices, mappings, min_count=None):
    """
    Utility function to find the number of matches between comparisons

    To check compatibility of data and mappings before running the entire "pipeline"
    Uses report matches function to print match count for all comparisons

    Args:
        nested_tags (dict): nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        int_matrices (dict): contains interaction matrices
            keys: sample-level tag of interaction matrix
            values: pd dataframe with pairwise interaction scores
        mappings (dict): id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        min_count (int, optional): Defaults to 10.
            match numbers lower than this value will be reported
    """
    # get comparisons
    tags = []
    for rep_tags in nested_tags.values():
        tags += rep_tags
    comps = list(combinations(tags, r=2))

    report_matches(
        comps,
        int_matrices,
        nested_tags,
        mappings,
        min_count=min_count)

### Functions to apply each analysis step separately to a set of profiles ###


def validate_params(arg_dict):
    """
    validates main function arguments

    Args:
        arguments (dict): main function's locals
                        used to get all arguments

    Raises:
        specific eror related to any problems with input arguments
    """

    # check output location
    output_location = arg_dict['output_location']
    if not os.path.isdir(output_location):
        msg = f'No such directory: {output_location}'
        raise NotADirectoryError(msg)

    # check top hit parameter compatibility
    th_criterium = arg_dict['th_criterium']
    th_percent = arg_dict['th_percent']
    if th_criterium == "percent":
        if not th_percent:
            msg = ("th_percent cannot be None if"
                   " th_criterium is percent")
            raise ValueError(msg)
    elif th_criterium == "best":
        pass
    else:
        msg = (f"{th_criterium} is not a valid value"
               " for th_criterium parameter."
               " choose 'percent' or 'best'")
        raise ValueError(msg)

    # check that replicate-level tags are unique
    sample_tags = ut.get_sample_tags(arg_dict['nested_tags'])
    if len(sample_tags) > len(set(sample_tags)):
        msg = ("replicate ids are not unique."
               " replicate ids must be unique across collections.")
        raise ValueError(msg)

    # check cluster annotation parameter compatibility
    coll_tags = arg_dict['nested_tags'].keys()
    perf_annot = arg_dict['perf_cluster_annotation']
    reference_groups = arg_dict['reference_groups']
    reference_tag = arg_dict['reference_tag']

    if perf_annot:
        if reference_groups is None:
            msg = ("if perf_cluster_annotation is True,"
                   " reference_groups cannot be None")
            raise ValueError(msg)
        if reference_tag is None:
            msg = ("if perf_cluster_annotation is True,"
                   " reference_tag cannot be None")
            raise ValueError(msg)
        elif reference_tag not in coll_tags:
            msg = (f"provided reference tag: {reference_tag} "
                   "not present in collection-level tags"
                   " provided in int_matrices param")
            raise ValueError(msg)


def process_params(nested_tags, job_name, output_location):
    """
    preprocessing of main function input parameters

    Args:
        nested_tags (dict): nested tag structure for profiles,
            keys: collection-level tags
            values: sample-level tags
        job_name (str): name of this job, used in output dir name
        output_location (str): folder path, output dir will be created here

    Raises:
        NotADirectoryError: if target location does not exists
        FileExistsError: if output folder name already exists

    Returns (tuple):
        job_name (string): name of this job, used in output dir name
        out_folder (string): path of folder to contain results
    """
    # check if location exists, create output folder
    if not job_name:
        job_name = '_'.join(list(nested_tags.keys()))

    # ADD NAME OF THE TOOL TO OUTPUT FOLDER, once it has a name
    out_folder = os.path.join(output_location, job_name + '_results')
    try:
        os.mkdir(out_folder)
    except FileExistsError as e:
        msg = f'file/folder with that name already exists: {out_folder}'
        raise FileExistsError(msg)

    return job_name, out_folder


def between_scoring(nested_tags, int_matrices, mappings,
                    p=0.90, min_search_weight=0.99,
                    th_criterium='percent', th_percent=1,
                    processes=1, chunksize=1000):
    """
    computes top hit scores between all interaction matrices

    first performs rbo scoring for all pairs between samples,
    then determines pairs that are reciprocal top hits

    Args:
        nested_tags (dict): dict with nested tag structure for samples
            keys: collection-level tags
            values: sample-level tags
        int_matrices (dict): contains interaction matrices
            keys: sample-level tag of interaction matrix
            values: pd dataframe with pairwise interaction scores
        mappings (dict): id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        p (float, optional): Defaults to 0.90.
            "top heaviness" parameter of rbo scoring
            lower values for p result in increased top-heaviness
        min_search_weight (float, optional): Defaults to 0.99.
            determines search depth of ranked lists. will search
            to a depth that results computation of fraction of total
            possible score equal to min_search_weight
        th_criterium (str, optional): Defaults to 'percent'.
            criterium to determine top hit
            "percent" counts proteins if both are in each
            other's top n % percent of ranked interactor lists
        th_percent (int, optional): Defaults to 1.
            top percent to consider when using "percent" th_criterium
        processes (int, optional): Defaults to 1.
            number of cpu cores to use
        chunksize (int, optional): Defaults to 1000.
            parallelization parameter

    Returns:
        dict: structure: {(left_tag)(right_tag):top_hits}
            top hits (pd series): top hits between two samples
    """

    # determine search depth
    eprint('determining search depth for rbo calculation..')
    collection_sizes = [matrix.shape[0] for matrix in int_matrices.values()]
    min_size = min(collection_sizes)
    search_depth = ps.det_search_depth(p, min_search_weight, min_size)
    eprint('search depth to reach min weight of'
           f' {min_search_weight}: {search_depth}')

    # get comparisons
    tags = []
    for rep_tags in nested_tags.values():
        tags += rep_tags
    comps = list(combinations(tags, r=2))

    # report the number of matches between comparisons before computing rbos
    eprint('checking the number of matching ids between comparisons..')
    report_matches(comps, int_matrices, nested_tags, mappings)
    eprint()

    # score each comparison, using score_comparison
    between_top_hits = {}
    for i, (left, right) in enumerate(comps):
        eprint('\rcomputing scores for comparison '
               f'{i+1} of {len(comps)}: {left}:{right}..', end="")
        left_scores = int_matrices[left]
        right_scores = int_matrices[right]

        mapping = ut.get_comp_mapping(left, right, nested_tags, mappings)

        rec_top_hits = score_comparison(
            left_scores, right_scores,
            p=p, mapping=mapping, search_depth=search_depth,
            th_criterium=th_criterium, th_percent=th_percent,
            processes=processes, chunksize=chunksize)

        between_top_hits[(left, right)] = rec_top_hits

    return between_top_hits


def score_comparison(left_scores, right_scores, mapping=False,
                     p=0.9, search_depth=None,
                     th_criterium='percent', th_percent=1,
                     processes=1, chunksize=1000):
    """
    Determine RBO scores and reciprocal top hits between pair of int matrices


    Args:
        [left|right]_scores (df): symmetric interaction matrix
            values: within-sample pairwise interaction scores
        mappings (dict): id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        p (float, optional): Defaults to 0.90.
            "top heaviness" parameter of rbo scoring
            lower values for p result in increased top-heaviness
        search_depth (_type_, optional): Defaults to None.
            number of ranks to consider when computing RBO scores
            if None, considers complete ranked lists
        th_criterium (str, optional): Defaults to 'percent'.
            criterium to determine top hit
            "percent" counts proteins if both are in each
            other's top n % percent of ranked interactor lists
        th_percent (int, optional): Defaults to 1.
            top percent to consider when using "percent" th_criterium
        processes (int, optional): Defaults to 1.
            number of cpu cores to use
        chunksize (int, optional): Defaults to 1000.
            parallelization parameter

    Returns:
        pd.Series: reciprocal top hits in pd series format
            2-level multiindex with id pair, values are scores
    """
    # determine rbo scores
    rbo_scores = ps.pairwise_rbo_scoring(
        left_scores, right_scores, mapping=mapping, p_param=p,
        search_depth=search_depth, processes=processes,
        chunksize=chunksize)

    # determine reciprocal top hits
    rec_top_hits = rh.get_reciprocal_top_hits(
        rbo_scores, score_type='between', criterium=th_criterium,
        percent=th_percent)

    # return reciprocal top hits
    return rec_top_hits


def get_within_top_hits(int_matrices, th_criterium='percent',
                        th_percent=1):
    """
    determines top hits within each interaction matrix

    Args:
        int_matrices (dict): contains interaction matrices
            keys: sample-level tag of interaction matrix
            values: pd dataframe with pairwise interaction scores
        th_criterium (str, optional): Defaults to 'percent'.
            criterium to determine top hit
            "percent" counts proteins if both are in each
            other's top n % percent of ranked interactor lists
        th_percent (int, optional): Defaults to 1.
            top percent to consider when using "percent" th_criterium

    Returns:
        dict: {'sample_tag':top_hits} structure
            top hits: pd series, top hits within each profile
    """

    within_top_hits = {}

    for tag, scores in int_matrices.items():
        print(f'determining reciprocal top hits within {tag} ..')
        within_top_hits[tag] = rh.get_reciprocal_top_hits(
            scores, score_type='within', criterium=th_criterium,
            percent=th_percent)
    return within_top_hits


def mcl_clustering(within_top_hits, between_top_hits,
                   out_folder, include_within=False,
                   wbratio=1, mcl_inflation=2,
                   processes=1):
    """
    perform mcl on network from a set of samples

    Args:
        [within|between]_top_hits (pd series):
            dict with within/between top hits per profile/comparison
        out_folder (string): 
            filepath of output directory
        include_within (bool, optional):Defaults to False.
            whether to include within profile interaction scores
            in combined network used as input for MCL clustering
        wbratio (int, optional): Defaults to 1.
            ratio of within/between score averages. within (interaction)
            and between (rbo) scores are normalised based on their average
        mcl_inflation (int, optional): Defaults to 2.
            inflation parameter of the mcl clustering algorithm. determines
            granularity of clustering result.
        processes (int, optional): Defaults to 1.
            number of cpu cores to use

    Returns:
        [mcl|network]_outfn (string):
            filepath of saved cluster and network results
    """
    # create combined network as input for MCL
    network_outfn = os.path.join(out_folder, 'combined_network.tsv')
    eprint(f'creating a combined network, saving to: {network_outfn}')
    mcl.create_combined_network(
        list(between_top_hits.values()),
        network_outfn,
        include_within=include_within,
        within_scores=list(within_top_hits.values()),
        wbratio=wbratio
    )

    # perform mcl clustering on combined network
    mcl_outfn = os.path.join(out_folder, 'mcl_result.tsv')
    eprint(f'performing mcl clustering, saving result to: {mcl_outfn}')
    mcl.run_MCL(
        network_outfn,
        mcl_outfn,
        inflation=mcl_inflation,
        processes=processes
    )
    eprint(f'clustering complete!')

    return mcl_outfn, network_outfn


def process_mcl_result(
        mcl_outfn, nested_tags, network_outfn, mappings,
        report_threshold=0.5,
        filter_clusters=True,
        perf_cluster_annotation=False,
        reference_groups=None,
        reference_tag=None,
        annot_fraction_threshold=0.5,
        annot_filter_mem_threshold=0.25):
    """
    processes and annotates raw mcl results for interpretation

    separate clusters per collection, computes cluster metrics

    Args:
        mcl_outfn (string):
            filepath of mcl result
        nested_tags (dict of dicts):
            nested tag structure for profiles
            keys: collection-level tags
            values: sample-level tags
        network_outfn (string):
            filepath of combined network
        mappings (dict):
            containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        report_threshold (float, optional): 
            in reporting mcl results of multi-sample collections,
            proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold. Defaults to 0.5.
        filter_clusters (bool, optional): default True
            whether to filter out clusters with less than 2 matches
            or proteins over report_threshold
        perf_cluster_annotation (bool, optional): Defaults to False.
            whether to perform automatic annotation of clusters using ref
        reference_groups (dict, optional): Defaults to None.
            only used if perf_cluster_annotation == True.
            contains the names(keys) and members(values) of the reference groups
        reference_tag (string, optional): Defaults to None.
            only used if perf_cluster_annotation == True.
            tag of collection on which annotation is to be based.
            member names in collection should match member names in reference
        annot_fraction_threshold (float, optional): Defaults to 0.5.
            minimum fraction of reference that should be present in
            the cluster to get as assignment
        annot_filter_mem_threshold (float, optional): Defaults to 0.25.
            cluster members scoring below threshold will be ignored
            in determining overlap of cluster with reference.
            if value is  None no filtering is applied

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
    # process the mcl results
    mcl_res = mcl.process_annot_MCL_res(
        mcl_outfn, nested_tags, network_outfn, mappings,
        report_threshold=report_threshold,
        filter_clusters=filter_clusters)

    # optionally annotate the clusters
    if perf_cluster_annotation:
        eprint('annotating clusters using reference..')

        # get subclusters
        subclusters = mcl_res['clusts_split'][reference_tag]

        # get cluster assignments
        assignments = annotate_clusters(
            subclusters, reference_groups,
            fraction_threshold=annot_fraction_threshold,
            filter_mem_threshold=annot_filter_mem_threshold)

        # add annotations to clust info
        mcl_res['clust_info']['reference_assignment'] = mcl_res['clust_info'].index.map(
            assignments)

    return mcl_res


def save_results(mcl_res, out_folder, mappings):
    """
    process and write human-readable results to file

    Args:
        mcl_res (dict): 
            dictionary with processed clustering results
        out_folder (string):
            path of output directory
        mappings (dict):
            containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
    """
    # clust info
    clust_info_outfn = os.path.join(out_folder, 'clust_info.tsv')
    mcl_res['clust_info'].to_csv(clust_info_outfn, sep='\t')

    # nodes
    nodes_outfn = os.path.join(out_folder, 'clust_nodes.tsv')
    mcl_res['nodes'].to_csv(nodes_outfn, sep='\t')

    # edges
    edges_outfn = os.path.join(out_folder, 'clust_edges.tsv')
    mcl_res['edges'].to_csv(edges_outfn, sep='\t')

    # get cluster member tables per collection
    member_tables = prd.split_clustmember_tables(mcl_res['nodes'], mappings)
    for name, table in member_tables.items():
        best_guesses = mcl_res['best_guess'][name]
        matches_over_threshold = mcl_res['match_over_threshold'][name]
        table['best_guess_selection'] = False
        table['match_over_threshold'] = False

        # add best guess selection to clustmember tables
        for cid, members in best_guesses.items():
            if cid not in table.index:
                continue
            selected = (table.index == cid) & (table['id'].isin(members))
            table.loc[selected, 'best_guess_selection'] = True

        # add matches over threshold to clustmember tables
        for cid, members in matches_over_threshold.items():
            if cid not in table.index:
                continue
            selected = (table.index == cid) & (table['id'].isin(members))
            table.loc[selected, 'match_over_threshold'] = True

        # write table to file
        table_outfn = os.path.join(out_folder, f'{name}_cluster_members.tsv')
        table.to_csv(table_outfn, sep='\t')

### Main function to perform complete analysis in 1 go ###


def main(nested_tags, int_matrices, mappings, p=0.90, min_search_weight=0.99,
         th_criterium='percent', th_percent=1,
         include_within=False,
         wbratio=1,
         mcl_inflation=2,
         output_location='.',
         job_name=None,
         report_threshold=0.5,
         filter_clusters=True,
         save_rthits=False,
         perf_cluster_annotation=False,
         reference_groups=None,
         reference_tag=None,
         annot_fraction_threshold=0.5,
         annot_filter_mem_threshold=0.25,
         processes=1, chunksize=1000):
    """
    complete rbo and clustering analysis from interaction matrices

    Args:
        nested_tags (dict): dict with nested tag structure for samples
            keys: collection-level tags
            values: sample-level tags
        int_matrices (dict): contains interaction matrices
            keys: sample-level tag of interaction matrix
            values: pd dataframe with pairwise interaction scores
        mappings (dict): id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        p (float, optional): Defaults to 0.90.
            "top heaviness" parameter of rbo scoring
            lower values for p result in increased top-heaviness
        min_search_weight (float, optional): Defaults to 0.99.
            determines search depth of ranked lists. will search
            to a depth that results computation of fraction of total
            possible score equal to min_search_weight
        th_criterium (str, optional): Defaults to 'percent'.
            criterium to determine top hit
            "percent" counts proteins if both are in each
            other's top n % percent of ranked interactor lists
            "best" takes only the single best hit
        th_percent (int, optional): Defaults to 1.
            top percent to consider when using "percent" th_criterium
        include_within (bool, optional): Defaults to False.
            whether to include within sample interaction scores
            in combined network used as input for MCL clustering
        wbratio (int, optional): Defaults to 1.
            ratio of within/between score averages. within (interaction)
            and between (rbo) scores are normalised based on their average
        mcl_inflation (int, optional): Defaults to 2.
            inflation parameter of the mcl clustering algorithm. determines
            granularity of clustering result.
        output_location (str, optional): Defaults to current working directory.
                filepath of output dir to be created
        job_name: str, default: concatenated collection tags
            name of this job, used in output dir name
        save_rthits (bool, optional): Defaults to False.
            whether to save reciprocal top hits to disk
        report_threshold (float, optional): in reporting mcl results of
            sampled collections, proteins are counted as member if they are
            a member of the cluster in a fraction of the
            samples >= report_threshold. Defaults to 0.5.
        filter_clusters (bool, optional): default True
            whether to filter out clusters with less than 2 matches
            or proteins over report_threshold
        perf_cluster_annotation (bool, optional): Defaults to False.
            whether to perform automatic annotation of clusters using ref
        reference_groups (dict, optional): Defaults to None.
            only used if perf_cluster_annotation == True.
            contains the names(keys) and members(values) of the reference groups
        reference_tag (string, optional): Defaults to None.
            only used if perf_cluster_annotation == True.
            tag of collection on which annotation is to be based.
            member names in collection should match member names in reference
        annot_fraction_threshold (float, optional): Defaults to 0.5.
            minimum fraction of reference that should be present in
            the cluster to get as assignment
        annot_filter_mem_threshold (float, optional): Defaults to 0.25.
            cluster members scoring below threshold will be ignored
            in determining overlap of cluster with reference.
            if value is  None no filtering is applied
        processes (int, optional): Defaults to 1.
            number of cpu cores to use

    Returns (dict): mcl_res, contains:
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
    # check if mcl is available, error if not
    if not ut.mcl_available():
        msg = ("MCL not available, cannot perform analysis"
               " if MCL executable is not in PATH")
        raise ValueError(msg)

    ## validate input parameters ##
    validate_params(locals())

    ## process input parameters ##
    job_name, out_folder = process_params(nested_tags, job_name,
                                          output_location)

    ## rename duplicates in the interaction data ##
    for name in int_matrices.keys():
        n_dups = prd.rename_duplicates_int_matrix(int_matrices[name])
        if n_dups > 0:
            print(f'{n_dups} duplicate ids found in {name},'
                  ' renamed occurences after first')

    ## between scoring ##
    between_top_hits = between_scoring(
        nested_tags,
        int_matrices,
        mappings,
        p=p,
        min_search_weight=min_search_weight,
        th_criterium=th_criterium,
        th_percent=th_percent,
        processes=processes,
        chunksize=chunksize)

    # add sample tags
    for (left, right), hits in between_top_hits.items():
        mindex = hits.index
        new_mindex = prd.add_tag_multiindex(left, right, mindex)
        between_top_hits[(left, right)].index = new_mindex

    # within top hits
    if include_within:
        within_top_hits = get_within_top_hits(
            int_matrices, th_criterium=th_criterium,
            th_percent=th_percent)

        # add sample tags
        for tag, hits in within_top_hits.items():
            mindex = hits.index
            new_mindex = prd.add_tag_multiindex(tag, tag, mindex)
            within_top_hits[tag].index = new_mindex

    else:
        within_top_hits = {}

    # optionally write top hits to file
    if save_rthits:
        eprint('writing reciprocal top hits to file')
        for tag, top_hit in within_top_hits.items():
            fn = os.path.join(out_folder, f'{tag}_within_top_hits.tsv')
            rh.save_top_hits(top_hit, fn)

        for cur_tags, top_hit in between_top_hits.items():
            fn = os.path.join(
                out_folder,
                f'{cur_tags[0]}|{cur_tags[1]}_between_top_hits.tsv')
            rh.save_top_hits(top_hit, fn)

    ## perform mcl clustering on the profile's top hits ##
    mcl_outfn, network_outfn = mcl_clustering(
        within_top_hits, between_top_hits,
        out_folder, include_within=include_within,
        wbratio=wbratio, mcl_inflation=mcl_inflation,
        processes=processes)

    eprint('processing clustering results..')
    mcl_res = process_mcl_result(
        mcl_outfn, nested_tags, network_outfn,
        mappings, report_threshold=report_threshold,
        filter_clusters=filter_clusters,
        perf_cluster_annotation=perf_cluster_annotation,
        reference_groups=reference_groups,
        reference_tag=reference_tag,
        annot_fraction_threshold=annot_fraction_threshold,
        annot_filter_mem_threshold=annot_filter_mem_threshold,
    )

    # write (human-readable) results to file
    save_results(mcl_res, out_folder, mappings)

    eprint(f'analysis complete! your results can be found here: {out_folder}')

    return mcl_res
