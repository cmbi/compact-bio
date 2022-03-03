try:
    import requests
except:
    msg = "cannot import requests package, download_sample_abuns function not available"
    print(msg)

def download_sample_abuns(
        sample_id,out_fn,
        output_ids = ['prot_ids'],
        prot_ids = [],
        id_type = 'prot_ids',
        url='https://www3.cmbi.umcn.nl/cedar/api/abundances'):
    """
    fetch abundances of a sample from cedar

    Args:
        sample_id (int): CEDAR CRS number of a sample
        out_fn (str): filepath, output location
        output_ids (list of strings, optional): Defaults to ['prot_ids'].
            the types of protein ids to include in the output.
            options: "prot_ids","prot_names","gene_names"
        prot_ids (list, optional): _description_. Defaults to [].
            if empty: fetches complete complexome profile
            otherwise: only fetch abundances for proteins matching given ids
        id_type (str, optional): Defaults to 'prot_ids'.
            type of identifier used when providing prot_ids
        url (str, optional): Defaults to 'https://www3.cmbi.umcn.nl/cedar/api/abundances'.
            url of CEDAR fetch_abundances api endpoint
    """
    headers={
        'Content-Type':'application/json',
        'Accept':'text/csv',
    }
    data = {
        'sample_id':sample_id,
        'prot_ids':prot_ids,
        'id_type':id_type,
        'output_ids':output_ids,
    }
    response = requests.post(url,json=data,headers=headers)
    if response.ok:
        print(f'response ok, writing result to: {out_fn} ')
        with open(out_fn, 'wb') as f_obj:
            f_obj.write(response.content)

def map_df_index(df,mapping):
    """
    rename df's index using given mapping {index:new_id}

    for ids that have no mapping original id is used
    """
    df = df.copy()
    df['mapped'] = df.index.map(mapping)

    # take care of missing labels
    missing = df['mapped'].isna()
    df.loc[missing,'mapped'] = df.loc[missing].index.values

    df.set_index('mapped',inplace=True)
    return df

def get_stripped_mapping(full_id_list,sep='::'):
    """
    get dict with stripped ids mapping to list of their full ids
    """
    mapping = {}

    for full_id in full_id_list:
        stripped = full_id.rsplit(sep,1)[0]
        if stripped in mapping.keys():
            mapping[stripped].append(full_id)
        else:
            mapping[stripped] = [full_id]
    return mapping

def get_cluster_max_fraction(clusters, profile):
    """
    determine fraction of cluster where mean at max
    """
    # scale profile to weigh each protein equally
    scaled = profile.scale()

    maxfracs = {}
    for name,members in clusters.items():
        frac = scaled.loc[members].mean().reset_index(drop=True).idxmax()
        maxfracs[name] = frac

    return maxfracs