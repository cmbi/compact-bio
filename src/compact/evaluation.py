#!/usr/bin/env
"""
metrics for evaluating clustering results

ari_ami: adjusted rand index and adjusted mutual information
        - to evaluate similarity between partitionings of the same protein set
mmr: maximum matching ratio
    - to evaluate clustering result against gold standard reference
geometric_accuracy: geometric mean of sensitivity and positive predictive value
    - to evaluate clustering result against gold standard reference

other:
-sensitivity
-positive predictive value
-overlap (metric between two sets)
"""
# import statements
import sys
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import adjusted_mutual_info_score
import pandas as pd
import numpy as np
import operator

# functions

#----- adjusted rand index, adjusted mutual information -----#


def ari_ami(left, right, report_type='ARI'):
    """
    calculate adujsted rand index for 2 partitions of same starting set

    requires identical source gene sets to be clustered
    Args:
        left,right: instances of dict with identical genes
        type('ARI' or 'AMI'): whether to compute ARI or AMI
    Returns:
        float: adjusted Rand index or adjusted mutual information
    """
    l_labels = left.as_labels
    r_labels = right.as_labels
    if len(l_labels) != len(r_labels):
        raise ValueError('partitioned sets are not equal in size!')
    labels = pd.concat([l_labels, r_labels], axis=1, ignore_index=False)
    if len(labels) != len(l_labels):
        raise ValueError('gene names in partitioned sets do not match!')
    if report_type == 'ARI':
        result = adjusted_rand_score(labels[0], labels[1])
    elif report_type == 'AMI':
        result = adjusted_mutual_info_score(labels[0], labels[1])
    else:
        raise ValueError(
            f"type {report_type} not supported. only 'ARI' or 'AMI' allowed")
    return result

#----- k-cliques -----#
# MAYBE to be implemented

#------- MMR -------#
# calculate overlap of each complex against reference, for each reference


def mmr(predicted, reference, report_best_matches=False):
    """
    calculate maximum matching ratio for predicted and reference dict

    described in Nepusz et al., 2012

    Args:
        predicted,reference: dicts with
            keys: cluster id
            values: list of cluster members
         report_best_matches:
            True to return which cluster matches each reference cluster best
            dict with key: reference cluster id, value: best match cluster id
    Returns:
        float: maximum matching ratio value
    """
    if not predicted or not reference:
        return 0
    best_matches = {}
    max_overlap_scores = []
    for ref_name, ref_set in reference.items():
        overlap_scores = {
            name: overlap(
                ref_set,
                pred_set) for name,
            pred_set in predicted.items()}
        best_match = max(overlap_scores.items(), key=operator.itemgetter(1))
        best_matches[ref_name] = best_match
        max_overlap_scores.append(best_match[1])
    mmr = np.mean(max_overlap_scores)
    if report_best_matches:
        return mmr, best_matches
    elif report_best_matches == False:
        return mmr
    else:
        raise ValueError('report_best_matches must be either True or False!')


def overlap(set1, set2):
    """
    calculate overlap between two sets (described in Bader & Hogue 2003)

    Args:
        set1/set2 (set): sets to determine overlap between

    Returns:
        float: overlap between two given input sets
    """
    if len(set1) == 0 or len(set2) == 0:
        return 0
    intersection = set(set1) & set(set2)
    overlap = (len(intersection)**2) / (len(set1) * len(set2))
    return overlap

#----- geometric accuracy -----#


def get_count_table(predicted, reference):
    """
    get count table for occurence of reference proteins in predicted complexes

    used to calculate sensitivity and ppv metrics
    Args:
        predicted/reference (dict): predicted and reference clusters

    Returns:
        pd.DataFrame: table with occurence counts of reference in predicted
    """
    count_series = {}
    for ref_name, ref_set in reference.items():
        count_dict = {}
        for pred_name, pred_set in predicted.items():
            intersect = set(ref_set) & set(pred_set)
            count_dict[pred_name] = len(intersect)
        count_series[ref_name] = count_dict
    result = pd.DataFrame.from_dict(count_series, orient='index')
    # remove reference complexes for which no proteins were in predicted
    result = result[(result.T != 0).any()]
    return result


def sensitivity(predicted, reference, disregard_unclustered=True):
    """
    sensitivity based on Brohee and van Helden (2006)

    different from original:only considers proteins that were detected.
    does take into account ref proteins that were detected but unclustered

    Args:
        predicted/reference (dict): predicted and reference clusters
        disregard_unclustered (bool, optional): Defaults to True.
            wheter to include 'unclustered'  as cluster set in
            Sn computation. if True, still uses unclustered
            proteins to total number of reference proteins

    Returns:
        float: sensitivity metric
    """
    if not predicted or not reference:
        return 0
    count_table = get_count_table(predicted, reference)
    rowsum = count_table.sum(axis=1)
    if disregard_unclustered and 'unclustered' in count_table.columns:
        count_table.drop('unclustered', axis=1, inplace=True)
    sensitivities = count_table.divide(rowsum, axis=0)
    complex_wise = sensitivities.max(axis=1)
    # sum of reference complex members that are actually clustered
    clust_rowsum = count_table.sum(axis=1)
    total_clust = clust_rowsum.sum()
    if total_clust == 0:
        return 0
    else:
        weighted_mean = (complex_wise * clust_rowsum).sum() / total_clust
        return weighted_mean


def ppv(predicted, reference, disregard_unclustered=True):
    """
    Positive predictive value(PPV) as described in Brohee and van Helden (2006)

    Args:
        predicted/reference (dict): predicted and reference clusters
        disregard_unclustered (bool, optional): Defaults to True.
            wheter to include 'unclustered'  as cluster set in
            ppv computation. if True, still uses unclustered
            proteins to total number of reference proteins

    Returns:
        float: positive predictive value metric
    """
    if not predicted or not reference:
        return 0
    count_table = get_count_table(predicted, reference)
    if disregard_unclustered and 'unclustered' in count_table.columns:
        count_table.drop('unclustered', axis=1, inplace=True)
    colsum = count_table.sum()
    ppvs = count_table.divide(colsum)
    cluster_wise = ppvs.max()
    # sum of cluster members without unclustered
    total_clust = colsum.sum()
    if total_clust == 0:
        return 0
    else:
        weighted_mean = (cluster_wise * colsum).sum() / total_clust
        return weighted_mean


def geometric_accuracy(predicted, reference):
    """
    geometric accuracy, based on Brohee and van Helden(2006)

    Sensitivity is different from Brohee and van Helden(see sensitivity function)

    Args:
        predicted/reference (dict): predicted and reference clusters

    Returns:
        float: geometric accuracy
    """
    if not predicted or not reference:
        return 0
    sn = sensitivity(predicted, reference)
    pposv = ppv(predicted, reference)
    geom_acc = np.sqrt((sn * pposv))
    return geom_acc


if __name__ == "__main__":

    # load clusterings
    cA = {'A': [1, 2, 3, 4, 5], 'B': [6, 7, 8], 'C': [9, 10]}
    cB = {'a': [1, 2, 3, 4, 5], 'b': [6, 7, 8], 'c': [9, 10]}
    cC = {'a': [1, 2, 3, 4], 'b': [6, 7, 8], 'c': [5, 9, 10]}
    cD = {'a': [1, 6, 9], 'b': [2, 3, 7], 'c': [4, 5, 8, 10]}
    cE = {'a': [1, 2, 3, 4], 'b': [9, 10]}
    cF = {'a': [1, 2, 3, 4], 'b': [6, 7], 'c': [9, 10],
          'unclustered': [5, 8]}

    complexes = {
        'complex 1': [1, 2, 3, 4, 5, 6, 7],
        'complex 2': [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
        'complex 3': [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
        'complex 4': [39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50],
    }

    predictions = {
        'cluster 1': [1, 2, 3, 4, 5, 6, 7],
        'cluster 2': [8, 9, 10, 11, 12, 13],
        'cluster 3': [14, 15, 16, 17, 18, 19, 20, 21],
        'cluster 4': [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 39, 40, 41, 42],
        'cluster 5': [36, 37, 38, 43, 44, 45, 46, 47],
        'unclustered': [48, 49, 50],
    }

    print(geometric_accuracy(cF, cA))
    print(mmr(cF, cA))
