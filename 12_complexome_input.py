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
            'CRS86': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS86_genenames.tsv',
            'CRS50': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS50_genenames.tsv',
            'CRS48': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS48_genenames.tsv',
            'CRS25': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS25_genenames.tsv',
            'CRS24': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS24_genenames.tsv',
            'CRS23': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS23_genenames.tsv',
            'CRS22': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS22_genenames.tsv',
            'CRS17': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS17_genenames.tsv',
        },
        'PF3_GAM':{
            'GAM1': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_GAM1_processed_renamed.tsv',
            'GAM2': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_GAM2_processed_renamed.tsv',
            'GAM3': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_GAM3_processed_renamed.tsv',
            'GAM4': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_GAM4_processed_renamed.tsv',            
        },
        'PF3_AS':{
            'AS1': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_AS1_processed_renamed.tsv',
            'AS2': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_AS2_processed_renamed.tsv',
            'AS3': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_AS3_processed_renamed.tsv',
            'AS4': '/home/joeri/Documents/Apicomplexa_project/data/profiles/Plasmo_AS4_processed_renamed.tsv',            
        },
        'PF3_SCH':{
            'CRS70': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS70_processed.tsv',
            'CRS71': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS71_processed.tsv',
            'CRS72': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS72_processed.tsv',
            'CRS73': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS73_processed.tsv',
            'CRS74': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS74_processed.tsv',
            'CRS75': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS75_processed.tsv',
        },
        'BER_SCH':{
            'CRS58': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS58_processed.tsv',
            'CRS59': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS59_processed.tsv',
            'CRS60': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS60_processed.tsv',
            'CRS61': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS61_processed.tsv',
            'CRS62': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS62_processed.tsv',
            'CRS63': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS63_processed.tsv',
        },
        'KNO_SCH':{
            'CRS64': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS64_processed.tsv',
            'CRS65': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS65_processed.tsv',
            'CRS66': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS66_processed.tsv',
            'CRS67': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS67_processed.tsv',
            'CRS68': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS68_processed.tsv',
            'CRS69': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS69_processed.tsv',
        },
        'TOX':{
            'TOX':'/home/joeri/Documents/Apicomplexa_project/data/profiles/Toxoplasma_comptab_processed.tsv',
        },
        'AT_LF':{
            'CRS100': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS100_abundances.tsv',
            'CRS101': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS101_abundances.tsv',
            'CRS102': '/home/joeri/Documents/Apicomplexa_project/data/profiles/CRS102_abundances.tsv',
        },
        'AT_SD':{
            'CRS103': '/home/joeri/Documents/Apicomplexa_project/data/profiles/seedlings1_abundances.tsv',
            'CRS104': '/home/joeri/Documents/Apicomplexa_project/data/profiles/seedlings2_abundances.tsv',
            'CRS105': '/home/joeri/Documents/Apicomplexa_project/data/profiles/seedlings3_abundances.tsv',
        },
        'ANOST':{
            'AN_DDM': '/home/joeri/Documents/Apicomplexa_project/data/profiles/anost_ddm.tsv',
            'AN_DIGI': '/home/joeri/Documents/Apicomplexa_project/data/profiles/anost_digi.tsv',
        },
        'BOVIN':{
            'BOV139': '/home/joeri/Documents/Apicomplexa_project/data/profiles/bovin_CRS139.tsv',
            'BOV140': '/home/joeri/Documents/Apicomplexa_project/data/profiles/bovin_CRS140.tsv',
            'BOV141': '/home/joeri/Documents/Apicomplexa_project/data/profiles/bovin_CRS141.tsv',
            'BOV142': '/home/joeri/Documents/Apicomplexa_project/data/profiles/bovin_CRS142.tsv',
            'BOV143': '/home/joeri/Documents/Apicomplexa_project/data/profiles/bovin_CRS143.tsv',
            'BOV144': '/home/joeri/Documents/Apicomplexa_project/data/profiles/bovin_CRS144.tsv',
        },
       'YARLI':{
           'YAR_1': '/home/joeri/Documents/Apicomplexa_project/data/profiles/yarli_r1.tsv',
           'YAR_2': '/home/joeri/Documents/Apicomplexa_project/data/profiles/yarli_r2.tsv',
            'YAR_3': '/home/joeri/Documents/Apicomplexa_project/data/profiles/yarli_r3.tsv',
            'YAR_4': '/home/joeri/Documents/Apicomplexa_project/data/profiles/yarli_r4.tsv',
        }
}

"""
mapping data files
------------------

filenames of identifier mapping files between
any pair of collections that have different 
identifiers
"""

mapping_data = {
    ('HUM','PF3_SCH'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/hum_pf_mapping.tsv",
    ('HUM','PF3_GAM'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/hum_pf_mapping.tsv",
    ('HUM','PF3_AS'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/hum_pf_mapping.tsv",
    ('BER_SCH','HUM'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/ber_hum_mapping.tsv",
    ('KNO_SCH','HUM'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/kno_hum_mapping.tsv",
    ('TOX','HUM'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/toxo_hum_mapping.tsv",
    ('BER_SCH','PF3_SCH'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/ber_pf_mapping.tsv",
    ('BER_SCH','PF3_GAM'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/ber_pf_mapping.tsv",
    ('BER_SCH','PF3_AS'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/ber_pf_mapping.tsv",
    ('PF3_SCH','KNO_SCH'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/pf_kno_mapping.tsv",
    ('PF3_GAM','KNO_SCH'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/pf_kno_mapping.tsv",
    ('PF3_AS','KNO_SCH'):"/home/joeri/Documents/Apicomplexa_project/data/mappings/pf_kno_mapping.tsv",
    ('TOX','PF3_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/toxo_pf_mapping.tsv',
    ('TOX','PF3_GAM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/toxo_pf_mapping.tsv',
    ('TOX','PF3_AS'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/toxo_pf_mapping.tsv',
    ('BER_SCH','KNO_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/ber_kno_mapping.tsv',
    ('BER_SCH','TOX'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/ber_tggt_mapping.tsv',
    ('KNO_SCH','TOX'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/kno_tggt_mapping.tsv',
    ('AT_LF','HUM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_hum.tsv',
    ('AT_SD','HUM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_hum.tsv',
    ('AT_LF','TOX'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_tggt.tsv',
    ('AT_SD','TOX'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_tggt.tsv',
    ('AT_LF','KNO_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_kno.tsv',
    ('AT_SD','KNO_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_kno.tsv',
    ('AT_LF','BER_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_ber.tsv',
    ('AT_SD','BER_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_ber.tsv',
    ('AT_LF','PF3_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_pf.tsv',
    ('AT_SD','PF3_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_pf.tsv',
    ('AT_LF','PF3_GAM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_pf.tsv',
    ('AT_SD','PF3_GAM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_pf.tsv',
    ('AT_LF','PF3_AS'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_pf.tsv',
    ('AT_SD','PF3_AS'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_pf.tsv',
    ('ANOST','AT_LF'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_arab_mapping.tsv',
    ('ANOST','AT_SD'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_arab_mapping.tsv',
    ('ANOST','HUM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_hum_mapping.tsv',
    ('ANOST','TOX'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_tggt_mapping.tsv',
    ('ANOST','KNO_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_kno_mapping.tsv',
    ('ANOST','BER_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_ber_mapping.tsv',
    ('ANOST','PF3_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_pf_mapping.tsv',
    ('ANOST','PF3_GAM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_pf_mapping.tsv',
    ('ANOST','PF3_AS'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_pf_mapping.tsv',
    ('ANOST','BOVIN'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_bovin_mapping.tsv',
    ('BER_SCH','BOVIN'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/ANKA_bovin_mapping.tsv',
    ('BER_SCH','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/ANKA_yarli_mapping.tsv',
    ('ANOST','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/anost_yarli_mapping.tsv',
    ('AT_LF','BOVIN'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_bovin_mapping.tsv',
    ('AT_SD','BOVIN'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_bovin_mapping.tsv',
    ('AT_LF','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_yarli_mapping.tsv',
    ('AT_SD','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/arab_yarli_mapping.tsv',
    ('BOVIN','KNO_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_PKNH_mapping.tsv',
    ('BOVIN','PF3_GAM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_plasmo_mapping.tsv',
    ('BOVIN','PF3_AS'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_plasmo_mapping.tsv',
    ('BOVIN','PF3_SCH'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_plasmo_mapping.tsv',
    ('BOVIN','TOX'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_tggt_mapping.tsv',
    ('BOVIN','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_yarli_mapping.tsv',
    ('KNO_SCH','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/PKNH_yarli_mapping.tsv',
    ('PF3_GAM','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/plasmo_yarli_mapping.tsv',
    ('PF3_AS','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/plasmo_yarli_mapping.tsv',
    ('PF3_SCH','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/plasmo_yarli_mapping.tsv',
    ('TOX','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/tggt_yarli_mapping.tsv',
    ('BOVIN','HUM'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/bovin_human_mapping.tsv',
    ('HUM','YARLI'):'/home/joeri/Documents/Apicomplexa_project/data/mappings/human_yarli_mapping.tsv',
}
