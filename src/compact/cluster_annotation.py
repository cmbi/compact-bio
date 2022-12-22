import operator

def annotate_clusters(clusters, reference, fraction_threshold=0.5,
                      filter_mem_threshold=0.25):
    """
    provide annotation to clusters using reference complexes

    bases annotation on subclusters of a collection of choice

    Args:
        clusters (dict): (sub)cluster members for one collection
                            key: clust_ids, values: list of members
        reference (dict): named reference groups
                            key: group name, values: list of members
        fraction_threshold (float, optional): Defaults to 0.5.
            minimum fraction of reference that should be present in
            the cluster to get as assignment
        filter_mem_threshold (float or None, optional): Defaults to 0.1
            cluster members scoring below threshold will be ignored
            in determining overlap of cluster with reference. 
            if value is  None no filtering is applied 

    Returns:
        dict: reference assignments for each cluster
            key: cluster id
            value: string with ref name(s) and fraction(s) present
    """

    # filter cluster members based on fraction present
    if filter_mem_threshold != None:
        clusters = filter_clust_members(
            clusters, threshold=filter_mem_threshold)

    # find best matches for reference sets
    matches = match_ref_to_clust(clusters, reference)

    # filter matches based on fraction_threshold
    matches = {key: val for key, val in matches.items()
               if val[1] > fraction_threshold}

    # aggregate assignments per cluster
    assignments = {}

    for name, (clust_id, score) in matches.items():
        if clust_id in assignments:
            assignments[clust_id] += f',{name}({score:.2f})'
        else:
            assignments[clust_id] = f'{name}({score:.2f})'

    return assignments


def match_ref_to_clust(predicted, reference, skip_singles=True):
    """
    find best cluster matches for reference complexes

    (Kw)Args:
        predicted,reference: dictionaries
        skip_singles:bool,default True
            skip reference complexes with only 1 subunit
    Returns:
        maximum matching ratio value
        optionally also returns dict with best matches
    """
    if not predicted or not reference:
        return 0
    best_matches = {}
    for ref_name, ref_set in reference.items():
        if skip_singles and len(set(ref_set)) == 1:
            continue
        overlap_scores = {name: fraction_present(
            ref_set, pred_set) for name, pred_set in predicted.items()}
        best_match = max(overlap_scores.items(), key=operator.itemgetter(1))
        best_matches[ref_name] = best_match
    return best_matches


def fraction_present(reference_set, cluster_set):
    """
    compute fraction of the reference complex present in the cluster

    Args:
        reference_set (list-like): reference group members
        cluster_set (list-like): cluster group members

    Returns:
        float: fraction of reference members present in cluster
    """
    intersection = set(reference_set) & set(cluster_set)
    if len(intersection) == 0:
        fraction = 0
    else:
        fraction = len(intersection) / len(reference_set)
    return fraction


def filter_clust_members(clusters, threshold=0.25):
    """
    Removes cluster members with a score below threshold

    Args:
        clusters (dict): clusters with members and associated scores
        threshold (float, optional): min score for inclusion. Defaults to 0.2.

    Returns:
        dict: clusters with members below threshold filtered out
    """
    filtered_clusters = {}
    for name, members in clusters.items():
        cur = {key: val for key, val in members.items() if val >= threshold}
        filtered_clusters[name] = cur
    return filtered_clusters
