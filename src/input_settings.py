"""
#################################################
file to provide input settings for run_compact.py
#################################################
"""

"""
sample data files
-----------------

gives names and the locations of the input sample files
as well as the names of the collections they belong to
"""
sample_data = {
    'HUM':{
        # 'CRS86': '../data/CRS86_genenames.tsv',
        'CRS50': '../data/CRS50_genenames.tsv',

    },
    'PF3_GAM':{
        'GAM1': '../data/Plasmo_GAM1_processed_renamed.tsv',
        'GAM2': '../data/Plasmo_GAM2_processed_renamed.tsv',
        
    },
    # 'TOX':{
    #     'TOX':'../data/Toxoplasma_comptab_processed.tsv',
    # }
}

"""
mapping data files
------------------

filenames of identifier mapping files between
any pair of collections that have different 
identifiers  
"""

mapping_data = {
    ("HUM","PF3_GAM"):"../data/mappings/hum_pf_mapping.tsv",
    ("TOX","HUM"):"../data/mappings/toxo_hum_mapping.tsv",
    ("TOX","PF3_GAM"):"../data/mappings/toxo_pf_mapping.tsv",
}