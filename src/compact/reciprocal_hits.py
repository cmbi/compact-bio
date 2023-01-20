"""
module to find protein pairs that are reciprocal "top" hits
with options for different criteria for "top" hit cut-off

main functions:
    - get_reciprocal_top_hits: determine reciprocal top hits from rbo score matrix between two datasets
"""

def get_top_hits(scores, criterium, percent=None, omit_self=False):
    """
    given a indexed series of scores, returns top hits

    Args:
        q_id (str): 
            query protein identifier
        scores (pd series): 
            scores with ids for query protein
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
    q_id = scores.name
    if omit_self:
        scores.drop(q_id, inplace=True)
    if criterium == 'best':
        top_hits = [scores.sort_values(ascending=False).index[0]]
    elif criterium == 'percent':
        if not (isinstance(percent, int) or isinstance(percent, float)):
            msg = 'when criterium is "percent", provide numeric percent argument'
            raise ValueError(msg)
        n = round(len(scores) * (percent / 100))
        top_hits = set(scores.sort_values(ascending=False).index[:n])
    
    else:
        msg = f'criterium not implemented!: {criterium}'
        raise NotImplementedError(msg)
    
    as_bools = scores.index.isin(top_hits)
    
    return as_bools

def get_reciprocal_top_hits(scores,score_type='between',criterium='percent',
                                 percent=1,out_type='series'):
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
            if 'dict', output is dict.
            if anything else, output is a dict

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

    index_top_hits = scores.apply(
        get_top_hits,axis=1,args=[criterium],percent=percent,
            omit_self=omit_self,result_type='broadcast').astype(bool)
    column_top_hits = scores.apply(
        get_top_hits,axis=0,args=[criterium],percent=percent,
                omit_self=omit_self,result_type='broadcast').astype(bool)

    reciprocal = index_top_hits & column_top_hits
    reciprocal_top_hits = scores[reciprocal].stack().sort_values(ascending=False)

    if out_type == 'dict':
        return reciprocal_top_hits.to_dict()
    
    else:
        return reciprocal_top_hits


def save_top_hits(top_hits, fn):
    """
    writes reciprocal top hits to tsv file

    Args:
        top_hits (pd series): reciprocal top hits
        fn (string): filepath of output
    """
    top_hits.to_csv(fn, sep='\t', header=False)
