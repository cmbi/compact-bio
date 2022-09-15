from . import utils as ut

def select_members(clusts_split,mappings,threshold=0.5):
    """
    Select best-guess complex members using matches and presence

    Args:
        clusts_split (_type_): _description_
        mappings (_type_): _description_
        threshold (float, optional): _description_. Defaults to 0.5.

    Returns: selected_members, match_filtered
        both are dict of dicts: contains dict of subclusters for each collection
            subclusters: list of members that passed filtering
    """
    thresh_filtered = threshold_filter(
        clusts_split,threshold=threshold)
    match_filtered = match_filter(
        clusts_split,thresh_filtered,mappings
    )
    selected_members = get_filter_union(
        thresh_filtered,match_filtered
    )

    return selected_members,match_filtered

def threshold_filter(clusts_split,threshold=0.5):
    """filter clusters based on frac_present threshold"""
    filtered = {}
    for tag,subclusters in clusts_split.items():
        tag_filtered = {}
        for cid,members in subclusters.items():
            thresh_members = set([key for key,val in members.items()
                              if val > threshold])
            tag_filtered[cid] = thresh_members
        filtered[tag] = tag_filtered
    return filtered

def match_filter(clusts_split,to_compare_clusts,mappings):
    """_summary_

    Args:
        clusts_split (_type_): _description_
        to_compare_clusts (_type_): _description_
        mappings (_type_): _description_

    Returns:
        _type_: _description_
    """
    filtered = {}
    # for every collection
    for tag,subclusters in clusts_split.items():
        # for each cluster
        tag_filtered = {}
        for cid,members in subclusters.items():
            # get corresponding to_compare cluster for each OTHER collection
            match_set = set([])
            for to_comp_tag,comp_subclusts in to_compare_clusts.items():
                if to_comp_tag == tag: # skip if current tag
                    continue
                # skip if this comp_tag does not have the current cluster 
                if not cid in comp_subclusts:
                    continue
                comp_members = comp_subclusts[cid]
                mapping = ut.get_col_mapping(tag,to_comp_tag,mappings)
                matches = ut.get_comparison_matches(
                    members,comp_members,mapping=mapping)
                # convert to list to get keys if its adict
                matches = list(matches) 
                # add matches with this comp_col to match_set
                match_set.update(matches)
            tag_filtered[cid] = match_set
        filtered[tag] = tag_filtered
    return filtered

def get_filter_union(left,right):
    """_summary_

    Args:
        left (_type_): _description_
        right (_type_): _description_

    Returns:
        _type_: _description_
    """
    union = {}
    for tag, subclusters in left.items():
        tag_union = {}
        for cid,left_members in subclusters.items():
            right_members = right[tag][cid]
            combined = left_members | right_members
            if len(combined) > 0:
                tag_union[cid] =  combined
        union[tag] = tag_union
    return union