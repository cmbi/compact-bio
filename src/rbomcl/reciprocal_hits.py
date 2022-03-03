"""
find protein pairs that are reciprocal "top" hits
option for different criteria for "top" hit cut-off
"""

# base library  imports
from sys import argv

# third party library imports
import pandas as pd


def get_top_hits(q_id, scores, criterium, percent=None, omit_self=False):
    """
    given a indexed series of scores, returns top hits

    Args:
        q_id (str): query protein identifier
        scores (pd series): scores with ids for query protein
        criterium ('best' or 'percent'): top hit criterium type
            if 'best': returns only the single best hit
            if 'percent': takes the top n % of highest scoring hits
        percent (numeric, optional): top percentage Defaults to None.
            if criterium is percent, returns top n percent of hits
        omit_self (bool, optional): Defaults to False.
            if True, drops query protein from scores list
    Raises:
        ValueError: raised if criterium is "percent" but no
            valid percent parameter is provided
        NotImplementedError: if invalid criterium parameter
            is provided

    Returns:
        list of strings: identifiers of top hits for query protein
    """
    if omit_self:
        scores.drop(q_id, inplace=True)
    if criterium == 'best':
        top_hits = [scores.sort_values(ascending=False).index[0]]
    elif criterium == 'percent':
        if not (isinstance(percent, int) or isinstance(percent, float)):
            msg = 'when criterium is "percent", provide numeric percent argument'
            raise ValueError(msg)
        n = round(len(scores) * (percent / 100))
        top_hits = list(scores.sort_values(ascending=False).index[:n])

    else:
        msg = f'criterium not implemented!: {criterium}'
        raise NotImplementedError(msg)
    return top_hits


def top_hits_to_series(top_hits_dict):
    """
    converts dict with top hits and scores to pd series

    Args:
        top_hits_dict (dict): dictionary with reciprocal
            top hits. strucure: {('l_id','r_id'):score}

    Returns:
        pd series: reciprocal top hits in pd series format
            2-level multiindex with id pair, values are scores
    """
    as_tuples = [(left, right, score) for (left, right), score in
                 top_hits_dict.items()]
    left, right, score = zip(*as_tuples)
    mindex = pd.MultiIndex.from_tuples(zip(left, right))
    top_hits = pd.Series(
        score, index=mindex
    )
    top_hits = top_hits.sort_values(ascending=False)
    return top_hits


def get_reciprocal_top_hits(scores, score_type='between', criterium='percent',
                            percent=1, out_type='series'):
    """
    determine reciprocal top hits for given score matrix

    Args:
        scores (pd df): matrix with interaction scores
            in case of 'within' scores: symmetric matrix with
            same ids in index and columns.
        score_type (str, optional): Defaults to 'between'.
            type of interaction scores. either 'within' or 'between'
            'within' for interaction scores within a single sample
            'between' for interaction scores between two samples
        criterium ('best' or 'percent'): top hit criterium type
            if 'best': returns only the single best hit
            if 'percent': takes the top n % of highest scoring hits
        percent (numeric, optional): top percentage Defaults to 1.
            if criterium is percent, returns top n percent of hits
        out_type (str, optional): Defaults to 'series'.
            if 'series', output is pd series.
            if anything else output is a dict

    Raises:
        ValueError: when score_type parameter is invalid

    Returns:
        pd series: reciprocal top hits in pd series format
            2-level multiindex with id pair, values are scores
        OR
        top_hits_dict (dict): dictionary with reciprocal
            top hits. strucure: {('l_id','r_id'):score}

    """
    # omit query id from potential hits
    #  if determining top hits "within" a sample
    if score_type == "within":
        omit_self = True
    elif score_type == "between":
        omit_self = False
    else:
        msg = f'"{score_type}" score type invalid, choose "within" or "between"'
        raise ValueError(msg)

    reciprocal_top_hits = {}
    # loop over vertical index
    for q_id, q_scores in scores.iterrows():

        # determine query top hits
        q_top_hits = get_top_hits(q_id, q_scores, criterium=criterium,
                                  percent=percent, omit_self=omit_self)

        # determine if reciprocal top hits
        for h_id in q_top_hits:
            h_scores = scores.loc[:, h_id]

            h_top_hits = get_top_hits(h_id, h_scores, criterium=criterium,
                                      percent=percent, omit_self=omit_self)

            # if reciprocal top hit, add to results
            if q_id in h_top_hits:
                if score_type == 'between':
                    score = scores.loc[q_id, h_id]
                    reciprocal_top_hits[(q_id, h_id)] = score
                else:
                    keys = reciprocal_top_hits.keys()
                    f_pair_present = (q_id, h_id) in keys
                    r_pair_present = (h_id, q_id) in keys
                    if not f_pair_present and not r_pair_present:
                        score = scores.loc[q_id, h_id]
                        reciprocal_top_hits[(q_id, h_id)] = score

    if out_type == 'series':
        reciprocal_top_hits = top_hits_to_series(reciprocal_top_hits)

    return reciprocal_top_hits


def save_top_hits(top_hits, fn):
    """
    writes reciprocal top hits to tsv file

    Args:
        top_hits (pd series): reciprocal top hits
        fn (string): filepath of output
    """
    top_hits.to_csv(fn, sep='\t', header=False)


if __name__ == "__main__":

    score_fn = argv[1]
    res = reciprocal_top_hits(score_fn, criterium='percent', percent=0.5)
    print(len(res))
