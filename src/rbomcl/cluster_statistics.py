from itertools import combinations
import multiprocessing as mp
from re import L

import numpy as np
import pandas as pd

import utils as ut
import process_data as prd

"""
module to perform statistics on clusters

to be used on raw clusters, before aggregating replicates etc.


given a cluster of a certain composition
    - determine number of matches within this cluster
        - code to count total matches
        - requires mappings, nested tags
        - take inspiration from process_cluster code that does something similar
    - sample a number of clusters randomly with the same distribution over samples
    - determine null distribution of number of matches
    - compute z-score, p-value
        - what fraction of nulls >= the scores? (p-value)
            - this doesn't make any assumptions about the distribution of n_matches
    - Then you need to do multiple testing correction to get FDR
        - the code to convert p-vals into q-vals is already complete in ../Apicomplexa_project/code/significance.py. Take from there!
            - CREATE SEPARATE MODULE FOR (MP) FDR CORRECTION. CAN BE USEFUL IN MANY CASES!!!!

SO. STEPS:
    - Determine what exact input you need, in what structure.
    - Then go from there
"""

def invert_nested_tags(nested_tags):
    """creates sample to collection mapping"""
    inverted = {}
    for ctag,mtags in nested_tags.items():
        for mtag in mtags:
            inverted[mtag] = ctag
    return inverted

def fetch_tag_members(members,tag):
    """
    for list of tagged members,get stripped ids with this tag

    Args:
        members (_type_): _description_
        tag (_type_): _description_

    Returns:
        _type_: _description_
    """
    tag_members = [mid[len(tag)+1:] for mid in members
                   if tag in mid]
    return tag_members

def split_clusters(clusts,sample_tags):
    """
    split clusters into members per sample

    Args:
        clusts (_type_): _description_
        sample_tags (_type_): _description_

    Returns:
        dict of dicts: _description_
    """
    nested_clusters = {}
    for cid, members in clusts.items():
        cur_clust = {}
        for tag in sample_tags:
            sample_members = fetch_tag_members(members,tag)
            if sample_members:
                cur_clust[tag] = sample_members
        nested_clusters[cid] = cur_clust
    return nested_clusters

def get_sample_tags(nested_tags):
    sample_tags = []
    for tags in nested_tags.values():
        sample_tags += tags
    return sample_tags

def count_comp_matches(nested_cluster,mappings,nested_tags):
    """
    """
    sample_tags = get_sample_tags(nested_tags)
    comps = list(combinations(sample_tags, r=2))
    
    comp_counts = {}
    for left,right in comps:
        # check if comp is present
        if not (left in nested_cluster.keys() and right in nested_cluster.keys()):
            continue

        left_members = prd.remove_appendices(nested_cluster[left])
        right_members = prd.remove_appendices(nested_cluster[right])
        comp_mapping = ut.get_comp_mapping(left,right,nested_tags,mappings)                
        matches = ut.get_comparison_matches(
            left_members,right_members,mapping=comp_mapping
        )
        comp_counts[(left,right)] = len(matches)
    return comp_counts

def sum_counts(comp_counts,nested_tags):
    """
    sum per collection and total count 
    
    in a way that loops over the total set of comparisons only once
    """
    summed_counts = {col:0 for col in nested_tags.keys()}
    summed_counts['total'] = 0
    sample_to_col = invert_nested_tags(nested_tags)
    for (left,right),count in comp_counts.items():
        left_col = sample_to_col[left]
        right_col = sample_to_col[right]

        # check if comparison is within collection
        if left_col == right_col:
            # if so add count to collection total
            summed_counts[left_col] += count

        # add count to total count
        summed_counts['total'] += count            
    
    return summed_counts

def get_match_counts(nested_clusters,mappings,nested_tags):
    """
    count matches for every cluster
    """
    match_counts = {}
    for cid,clust in nested_clusters.items():
        print(f'\rcounting matches for cluster: {cid}',end="")
        clust_matches = count_comp_matches(clust,mappings,nested_tags)
        summed_counts = sum_counts(clust_matches,nested_tags)
        match_counts[cid] = summed_counts
    print('\n')
    return match_counts

def sample_null(composition,sample_ids):
    """
    sample null clusters with same sample composition as given cluster
    """
    null_sample = {}
    for sample,size in composition.items():
        s_ids = sample_ids[sample]
        sampled = np.random.choice(s_ids,size,replace=False)
        null_sample[sample] = sampled
    return null_sample

def sample_null_counts(nested_cluster,sample_ids,nested_tags,mappings,
                    n=1000):
    """
    sample null distribution of counts for given cluster
    """
    composition = {key: len(val) for key,val in nested_cluster.items()}

    null_counts = []
    for i in range(n):
        print(f'\rtaking null sample {i+1} of {n}',end="")
        null_sample = sample_null(composition,sample_ids)
        matches = count_comp_matches(null_sample,mappings,nested_tags)
        summed_counts = sum_counts(matches,nested_tags)
        null_counts.append(summed_counts)
    print('\n')
    return null_counts

def summarize_null(null_counts,real_counts):
    """summarizes null counts and pval for real counts"""
    df = pd.DataFrame(null_counts)
    means = df.mean().to_dict()
    stds = df.std().to_dict()
    n = df.shape[0]

    pvals = {}
    for comp,count in real_counts.items():
        if count == 0:
            p = 1
        else:
            p = ((df[comp] > count).sum())/n
        pvals[comp] = p

    return means,stds,pvals

def score_cluster(nested_cluster,real_count,sample_ids,nested_tags,
                  mappings,n):
    """
    score given cluster by sampling null and comparing to real counts
    """
    nulls = sample_null_counts(nested_cluster,sample_ids,nested_tags,
                               mappings,n=n)
    return summarize_null(nulls,real_count)

def score_cluster_named(name,args):
    return name,score_cluster(*args)

def score_clusters(nested_clusters,real_counts,sample_ids,nested_tags,
                   mappings,n=1000,processes=1):
    """
    apply cluster scoring to given set of clusters
    """

    to_iter = (
            (cid,
            (
                clust,
                real_counts[cid],
                sample_ids,
                nested_tags,
                mappings,n
            )
        )
        for cid,clust in nested_clusters.items())
    pool = mp.Pool(processes)
    result = pool.starmap(score_cluster_named,to_iter)

    cluster_means = {}
    cluster_stds = {}
    cluster_pvals = {}
    for name,(means,stds,pvals) in result:
        cluster_means[name] = means
        cluster_stds[name] = stds
        cluster_pvals[name] = pvals

    return cluster_means,cluster_stds,cluster_pvals

    for cid,clust in nested_clusters.items():
        print(f'\rscoring cluster: {cid}',end="")
        means,stds,pvals = score_cluster(clust,real_counts[cid],sample_ids,
                                         nested_tags,mappings,n=n)
        cluster_means[cid] = means
        cluster_stds[cid] = stds
        cluster_pvals[cid] = pvals
    print('\n')
    return cluster_means,cluster_stds,cluster_pvals

if __name__ == "__main__":    # THINK ABOUT WHAT STRUCTURE I SHOULD HAVE THE SAMPLES IN

    import process_data as prd
    from run_compact import parse_settings, parse_mappings,get_nested_tags,parse_profiles,get_int_matrices
    import pandas as pd

    mcl_res_fn = '/home/joerivs/Documents/Apicomplexa_project/results/c12_run_Apr1_results/mcl_result.tsv'
    settings_fn = '12_complexome_input.py'
    clusts = prd.parse_MCL_result(mcl_res_fn)
    sample_data,mapping_data = parse_settings(settings_fn)
    mappings = parse_mappings(mapping_data)
    profiles = get_int_matrices(parse_profiles(sample_data))

    sample_ids = {name:profile.index.values 
                  for name,profile in profiles.items()}

    # print(sample_ids.keys())

    # decide later how to get these in actual code
    nested_tags = get_nested_tags(sample_data)
    sample_tags = []
    for tags in nested_tags.values():
        sample_tags += tags

    # separate clusters into subclusters per sample
    nested_clusters = split_clusters(clusts,sample_tags)

    #### TAKE A SUBSET OF CLUSTERS FOR TESTING
    # to_test = [5,6,7,10,20,100]
    # nested_clusters = {i:nested_clusters[i] for i in to_test}
    ##########################################

    # count total and within-subcluster matches
    real_counts = get_match_counts(nested_clusters,mappings,nested_tags)

    # as_df = pd.DataFrame.from_dict(res,orient='index')
    # as_df.to_csv('~/Documents/Apicomplexa_project/results/c12_run_Apr1_results/cluster_match_counts.tsv',sep='\t')

    #filter relevant clusters
    print('filtering clusters..')
    filtered_ids = [cid for cid,counts in real_counts.items()
                    if counts['total'] >= 2]
    print(len(real_counts))
    print(len(filtered_ids))
    filtered_clusters = {i:nested_clusters[i] for i in filtered_ids}
    filtered_counts = {i:real_counts[i] for i in filtered_ids}

    # compute null-based scores for relevant clusters
    print('scoring relevant clusters:')
    means,stds,pvals = score_clusters(filtered_clusters,filtered_counts,
                   sample_ids, nested_tags,mappings,n=1000,processes=5)

    means = pd.DataFrame.from_dict(means,orient='index')
    stds = pd.DataFrame.from_dict(stds,orient='index')
    pvals = pd.DataFrame.from_dict(pvals,orient='index')

    print(means.shape)
    print(stds.shape)
    print(pvals.shape)

    means.to_csv('/home/joerivs/Downloads/test_cluster_means.tsv',sep='\t')
    stds.to_csv('/home/joerivs/Downloads/test_cluster_stds.tsv',sep='\t')
    pvals.to_csv('/home/joerivs/Downloads/test_cluster_pvals.tsv',sep='\t')


    # print(as_df.shape)
    # print(as_df.head())

    # c0_matches = count_comp_matches(
    #     nested_clusters[10],
    #     mappings,
    #     nested_tags)


    # real_counts = sum_counts(c0_matches,nested_tags)

    # null_counts = sample_null_counts(
    #     nested_clusters[10],
    #     sample_ids,
    #     nested_tags,
    #     mappings,
    #     n=10)
    
    # print(null_counts)
    # null_df = pd.DataFrame(null_counts)
    # print(null_df)
    # print(null_df.describe())

    # print(real_counts)
    # print(real_counts)
    # print(len(null_counts))
    # print(res)
    # print(len(c0_matches))
    # df = pd.DataFrame.from_dict(c0_matches)
    # print(df.head(),orient='index')


    # print(pd.DataFrame.from_dict())
    # print(nested_clusters.keys())
    # print(nested_clusters[0].keys())