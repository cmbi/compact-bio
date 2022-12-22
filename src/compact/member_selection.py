from . import utils as ut


def select_members(clusts_split, mappings, threshold=0.5):
    """
    Select best-guess complex members using matches and presence

    Args:
        clusts_split (nested dict structure):
                cluster members, split per collection, samples aggregated
                members are dict with key: member id, value: member weight
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject
        threshold (float, optional): Defaults to 0.5.
            the fraction clustered threshold. members with a score
            GREATER than this threshold will be selected

    Returns: selected_members, match_filtered
        both are dict of dicts: contains dict of subclusters for each collection
            subclusters: list of members that passed filtering
    """
    thresh_filtered = threshold_filter(
        clusts_split, threshold=threshold)
    match_filtered = match_filter(
        clusts_split, thresh_filtered, mappings
    )
    selected_members = get_filter_union(
        thresh_filtered, match_filtered
    )

    return selected_members, match_filtered


def threshold_filter(clusts_split, threshold=0.5):
    """
    filter clusters based on frac_present threshold

    Args:
        clusts_split (dict):
                cluster members, split per collection, samples aggregated
        threshold (float, optional): Defaults to 0.5.
            the fraction clustered threshold. members with a score
            GREATER than this threshold will be selected

    Returns:
        dict of dict:
            clusts split, but only with members that pass threshold
            weights are not stored just identifiers
    """
    filtered = {}
    for tag, subclusters in clusts_split.items():
        tag_filtered = {}
        for cid, members in subclusters.items():
            thresh_members = set([key for key, val in members.items()
                                  if val > threshold])
            tag_filtered[cid] = thresh_members
        filtered[tag] = tag_filtered
    return filtered


def match_filter(clusts_split, to_compare_clusts, mappings):
    """
    filter to_compare_clusters based on match counts

    Args:
        clusts_split (dict):
                cluster members, split per collection, samples aggregated
        to_compare_clusts (dict of dicts):
            split clusters with member ids,
            but only with members that pass threshold
        mappings (dict): containing id mappings between collections
            keys: tuple with (query,subject) collection-level tags
            values: dicts with id mappings from query to subject

    Returns:
        dict of dict:
            clusts split filtered on members that have a match that
            passes threshold. weights are not stored just identifiers
    """
    filtered = {}
    # for every collection
    for tag, subclusters in clusts_split.items():
        # for each cluster
        tag_filtered = {}
        for cid, members in subclusters.items():
            # get corresponding to_compare cluster for each OTHER collection
            match_set = set([])
            for to_comp_tag, comp_subclusts in to_compare_clusts.items():
                if to_comp_tag == tag:  # skip if current tag
                    continue
                # skip if this comp_tag does not have the current cluster
                if cid not in comp_subclusts:
                    continue
                comp_members = comp_subclusts[cid]
                mapping = ut.get_col_mapping(tag, to_comp_tag, mappings)
                matches = ut.get_comparison_matches(
                    members, comp_members, mapping=mapping)
                # convert to list to get keys if its adict
                matches = list(matches)
                # add matches with this comp_col to match_set
                match_set.update(matches)
            tag_filtered[cid] = match_set
        filtered[tag] = tag_filtered
    return filtered


def get_filter_union(left, right):
    """
    given two nested cluster dicts, returns their union

    Args:
        left/right (dict of dicts):
            split clusters per collection with member id list

    Returns:
    dict of dicts: split clusters per collection with member id list
    """
    union = {}
    for tag, subclusters in left.items():
        tag_union = {}
        for cid, left_members in subclusters.items():
            right_members = right[tag][cid]
            combined = left_members | right_members
            if len(combined) > 0:
                tag_union[cid] = combined
        union[tag] = tag_union
    return union
