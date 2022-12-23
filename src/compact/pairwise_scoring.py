"""
module to perform pairwise rbo scoring between interaction datasets

Important functions:
    - pairwise_rbo_scoring: perform rbo scoring between two interaction matrices
    - score_comparison: rbo scoring between two correlation datasets
    - determine_top_weightedness: determine contribution of top d ranks to rbo score
    - det_search_deph: determine required search depth for given value for p parameter an mininum score weight

"""
# base imports
from itertools import product
import multiprocessing as mp

# third party library imports
import pandas as pd
import numpy as np

# local library imports
from . import process_data as prd
from .utils import get_comparison_matches
from .utils import eprint

# import available rbo implementation
try:
    from fastrbo import rank_biased_overlap

    def compute_rbo(search_depth, p, lranks, rranks):
        return rank_biased_overlap(
            search_depth, p, lranks, rranks)

except BaseException:
    eprint('fastrbo package not installed. using slower rbo package')
    eprint('install fastrbo to reduce computation time: https://github.com/joerivstrien/fastrbo')
    from rbo import RankingSimilarity

    def compute_rbo(search_depth, p, lranks, rranks):
        rbo_ranking = RankingSimilarity(lranks, rranks)
        return rbo_ranking.rbo(p=p)


def rename_indices(left, right, mapping):
    """
    assign new names to left and right indexes

    In a way that ensures that matching proteins
    will have matching names. left index should
    match with mapping keys, right index ids should
    match with mapping values.

    Args:
        left|right (list-like): indexes of compared samples
        mapping (dict):
            id mapping from left to right

    Returns:
        tuple (list,list): renamed left and right indexes
    """
    left_stripped = prd.remove_appendices(left)
    right_stripped = prd.remove_appendices(right)

    matches = get_comparison_matches(left_stripped, right_stripped,
                                     mapping=mapping)

    eprint(f'\nnumber of mapping ids found between profiles: {len(matches)}')

    renamed = assign_matched(left_stripped, right_stripped, matches)
    return renamed


def assign_matched(left, right, left_dict):
    """
    create renamed indices using matched ids

    Args:
        left|right (list-like): indexes of compared samples
        left_dict (dict): contains matching id pairs
            key: left id, value: right id

    Returns:
        tuple (list,list): renamed left and right index
    """
    # get inverse dictionary for renaming right index
    right_dict = prd.invert_mapping(left_dict)
    new_left = []
    for val in left:
        if val in left_dict.keys():
            new_val = val + '|' + left_dict[val]
            new_left.append(new_val)
        else:
            new_left.append(val)

    new_right = []
    for val in right:
        if val in right_dict.keys():
            new_val = right_dict[val] + '|' + val
            new_right.append(new_val)
        else:
            new_right.append(val)

    return new_left, new_right


def compute_ranked_lists(interaction_scores, search_depth=None):
    """
    compute ranked list for each protein based on interaction scores

    Args:
        interaction_scores (pd df): symmetric interaction matrix
            values: within-sample pairwise interaction scores
        search_depth (int or None, optional): Defaults to None.
            determines length of ranked list
            if None, all interactor ranks included

    Returns:
        list of lists: for each protein, contains a
        sorted ranked list of top interactors,
        not including self
    """
    ranked_lists = []
    for prot, scores in interaction_scores.iterrows():
        sorted_scores = scores.sort_values(ascending=False)
        sorted_list = list(sorted_scores.index.values)
        # remove (only first occurence) of current protein id
        # in case of duplicate id's assumes
        # current protein has higher score than duplicate
        sorted_list.remove(prot)

        # in case of duplicate ids in ranked list,
        # rename lowest scoring protein, so highest
        # gets considered by scoring algorithm
        sorted_list = prd.rename_duplicates(sorted_list)

        if search_depth:
            sorted_list = sorted_list[:search_depth]

        ranked_lists.append(sorted_list)

    return ranked_lists


def compute_rbo_scores(ranked_left, ranked_right, p,
                       processes=1,
                       chunksize=1000):
    """
    compute RBO scores for all pairs between samples

    Args:
        ranked_[left|right] (list of tuples):
            list of ranked lists for each protein in respective sample
            tuple structure: (prot_id,ranked_list)
        p (float): Rank Biased Overlap "p" parameter range: 0 to 1
            this parameter determines top-weightedness of rbo metric
            lower values result in more top-weightedness
        processes (int, optional): number of processes/threads. Defaults to 1.

    Returns:
        list of tuples: rbo scores for each protein pair
        tuple structure: ((left_id,right_id),rbo_score)
    """

    # create iterable that goes over all pairs
    comparisons = product(ranked_left, ranked_right)

    # compute rbo for each pair, using multiprocessing
    score_list = score_comparisons_parallel(comparisons, p,
                                            processes=processes,
                                            chunksize=chunksize)
    return score_list


def score_comparisons_parallel(comparisons, p,
                               processes=mp.cpu_count() - 1,
                               chunksize=1000):
    """
    wrapper function to paralellize score_comparison

    Args:
        comparisons (iterator): pairs of ranked lists to be compared
            ranked_list (tuple): (prot_id,ranked_list)
        p (float): Rank Biased Overlap "p" parameter range: 0 to 1
            this parameter determines top-weightedness of rbo metric
            lower values result in more top-weightedness
        processes (int, optional): number of processes/threads. Defaults to 1.

    Returns:
        list of tuples: rbo scores for each protein pair
        tuple structure: ((left_id,right_id),rbo_score)
    """
    pool = mp.Pool(processes)
    to_iter = ((comp, p) for comp in comparisons)
    result = pool.starmap(score_comparison, to_iter,
                          chunksize=chunksize)
    return result


def score_comparison(comparison, p):
    """
    compute rbo score between two given ranked lists

    Args:
        comparison (tuple): pair of ranked lists to be compared
            ranked_list (tuple): (prot_id,ranked_list)
        p (float): Rank Biased Overlap "p" parameter range: 0 to 1
            this parameter determines top-weightedness of rbo metric
            lower values result in more top-weightedness

    Returns:
        tuple: containing pair of ids and rbo score
            structure: ((left_id,right_id),rbo_score)
    """
    (lprot, lranks), (rprot, rranks) = comparison
    rbo = compute_rbo(len(lranks), p, lranks, rranks)
    return ((lprot, rprot), rbo)


def result_to_df(result, left_index, right_index,
                 structure='matrix'):
    """
    convert list-format rbo scoring result to table

    Args:
        result (list of tuples): rbo scores for each protein pair
            tuple structure: ((left_id,right_id),rbo_score)

        [left|right]_index (pd Index):
            index of left and right sample to reorder resulting dataframe
        structure ("matrix" or "sorted_list", optional): output format.
            Defaults to 'matrix'.
            if "matrix": matrix-like dataframe with:
                rows: left index. columns: right index, values are scores
            if "sorted_list": long-form list with score pairs
                columns: left_id,right_id,score
    Raises:
        ValueError: in case of invalid structure parameter

    Returns:
        pd dataframe: rbo scoring result in pd df format
            either in matrix or sorted list structure
    """
    if structure == 'matrix':
        pair, score = zip(*result)
        left, right = zip(*pair)
        df = pd.DataFrame({
            'left_id': left,
            'right_id': right,
            'score': score,
        })
        df.set_index(['left_id', 'right_id'], inplace=True)

        df = df.unstack()
        df.columns = df.columns.droplevel(0)
        df = df.loc[left_index, right_index]
    elif structure == 'sorted_list':
        sorted_result = sorted(result, key=lambda x: x[1],
                               reverse=True)
        pair, score = zip(*sorted_result)
        left, right = zip(*pair)
        df = pd.DataFrame({
            'left_id': left,
            'right_id': right,
            'score': score,
        })
    else:
        msg = f'"{structure}" is not a valid structure parameter'
        raise ValueError(msg)

    return df


def pairwise_rbo_scoring(left_scores, right_scores, mapping=False,
                         p_param=0.90, search_depth=None,
                         processes=1, chunksize=1000):
    """
    perform rbo scoring between 2 interaction matrices

    Args:
        [left|right]_scores (pd df): symmetric interaction matrix
            values: within-sample pairwise interaction scores
        mapping (bool, optional): Defaults to False.
            if not False: mapping (dict): id mapping/orthology between two collections
            keys: identifiers of "left" collection
            values: corresponding identifiers of "right" collection
        p_param (float, optional): Defaults to 0.90.
            Rank Biased Overlap "p" parameter range: 0 to 1
            this parameter determines top-weightedness of rbo metric
            lower values result in more top-weightedness
        search_depth (int, optional): Defaults to None.
            number of ranks to consider when computing RBO scores
            if None, considers complete ranked lists
        processes (int, optional): Defaults to 1.
            number of processes/threads. 
    Returns:
        pd df: rbo scores for all pairs between left right
            matrix-structured dataframe with:
            rows: left index. columns: right index, values are rbo scores
    """

    old_left_index = left_scores.index
    old_right_index = right_scores.index

    # create a copy of the score frames as to not alter
    # the index/columns of the dfs provided
    left_scores = left_scores.copy()
    right_scores = right_scores.copy()

    # rename protein indices, matching mapping ids
    # ASSUMING scores dfs have same order in index and columns!
    if mapping:
        new_left_index, new_right_index = rename_indices(
            left_scores.index, right_scores.index, mapping)
        left_scores.index = new_left_index
        right_scores.index = new_right_index
        left_scores.columns = new_left_index
        right_scores.columns = new_right_index
    else:
        # do not need to make ids match using mapping
        # but do need to remove appendages to not miss matches
        left_stripped = prd.remove_appendices(old_left_index.values)
        right_stripped = prd.remove_appendices(old_right_index.values)
        left_scores.index = left_stripped
        left_scores.columns = left_stripped
        right_scores.index = right_stripped
        right_scores.columns = right_stripped

    # get ranked lists for each protein in both interaction matrices
    # uses renamed indices where matched proteins have the same ids

    left_ranked_lists = compute_ranked_lists(
        left_scores, search_depth=search_depth)
    right_ranked_lists = compute_ranked_lists(
        right_scores, search_depth=search_depth)

    # assign prot ids to each ranked list by zipping in old index before
    # renaming
    left_ranked_lists = list(zip(old_left_index.values, left_ranked_lists))
    right_ranked_lists = list(zip(old_right_index.values, right_ranked_lists))

    # compute rbo scores for protein pairs between samples
    rbo_scores = compute_rbo_scores(left_ranked_lists,
                                    right_ranked_lists,
                                    p=p_param,
                                    processes=processes,
                                    chunksize=chunksize)

    # provide old indices, to sort the rbo result indices
    # to match the order in the input data
    score_df = result_to_df(rbo_scores, old_left_index,
                            old_right_index)

    return score_df


def determine_top_weightedness(p, d):
    """
    determine the contribution the top d ranks have to the rbo score

    taken from https://github.com/changyaochen/rbo/blob/master/rbo/rbo.py

    Args:
        p (float): Rank Biased Overlap "p" parameter range: 0 to 1
            this parameter determines top-weightedness of rbo metric
            lower values result in more top-weightedness
        d (int): 
            top d ranks for which to determine contribution

    Returns:
        float: the fraction of the rbo score that is determined
            by the top d ranks for the given p parameter
    """
    if d == 0:
        top_w = 1
    elif d == 1:
        top_w = 1 - 1 + 1.0 * (1 - p) / p * (np.log(1.0 / (1 - p)))
    else:
        sum_1 = 0
        for i in range(1, d):
            sum_1 += 1.0 * p**(i) / i
        top_w = 1 - p**(i) + 1.0 * (1 - p) / p * (i + 1) * \
            (np.log(1.0 / (1 - p)) - sum_1)  # here i == d-1
    return top_w


def det_search_depth(p, min_weight, shortest_list_len,
                     stepsize=1):
    """
    determine required search depth for given p and min_weight

    Args:
        p (float): Rank Biased Overlap "p" parameter range: 0 to 1
            this parameter determines top-weightedness of rbo metric
            lower values result in more top-weightedness
        min_weight (float): minimal weight of contribution to total
            rbo score requirement
        shortest_list_len (int): length of shortest ranked list to be
            used in calculating rbo score (search depth will be maximally
            this length)
        stepsize (int, optional): Defaults to 1.
            stepsize for walking over search depth values to consider.
            will stop if it reaches search depth with adequate weight contribution

    Returns:
        int: the selected search depth
    """
    # cannot compute weight if p =1, but will usually need entire list anyway
    if p == 1.0:
        return shortest_list_len
    min_achieved = False
    cur_rank = 1
    while min_achieved == False:
        if cur_rank >= shortest_list_len:
            eprint('min_weight not achievable,'
                   ' returning shortest_list_len')
            return shortest_list_len
        fraction = determine_top_weightedness(p, cur_rank)
        if fraction >= min_weight:
            min_achieved = True
        else:
            cur_rank += stepsize

    return cur_rank
