from compact import MCL_clustering as mcl
import pytest
import pandas as pd
from pandas.testing import assert_series_equal,assert_frame_equal

@pytest.fixture
def score_list():
    return [
        pd.Series({
            'prot1':0.9,
            'prot2':0.8,
            'prot3':0.7,
        }),
        pd.Series({
            'prot4':0.5,
            'prot5':0.3,
            'prot6':0.4,
        })
    ]

@pytest.fixture
def clusters():
    return {
        0:[
            'S1_prot1',
            'S2_prot1',
            'S1_prot3',
            'S3_prot1',
            'S4_prot1',
        ],
        1:[
            'S1_prot4',
            'S3_prot2',
            'S5_prot6',
            'S4_prot2',
        ],
    }

@pytest.fixture
def separated_clusters():
    return {
        'col1': {0: {'prot1': 1.0, 'prot3': 0.5}, 1: {'prot4': 0.5}},
        'col2': {0: {'prot1': 1.0}, 1: {'prot2': 1.0}},
        'col3': {1: {'prot6': 1.0}}
    }

@pytest.fixture
def sample_split_clusters():
        return {
        0:{
            'S1':['prot1','prot3'],
            'S2':['prot1'],
            'S3':['prot1'],
            'S4':['prot1'],
        },
        1:{
            'S1':['prot4'],
            'S3':['prot2'],
            'S4':['prot2'],
            'S5':['prot6'],
        }
    }


@pytest.fixture
def nested_tags():
    return {
        'col1':['S1','S2'],
        'col2':['S3','S4'],
        'col3':['S5'],
    }

@pytest.fixture
def inverted_tags():
    return {
        'S1':'col1',
        'S2':'col1',
        'S3':'col2',
        'S4':'col2',
        'S5':'col3',
    }

@pytest.fixture
def mappings():
    return {
        ('col1','col3'):{
            'prot4':'prot6'
        },
    }

@pytest.fixture
def counts():
    return pd.DataFrame({
        'col1': {0: 1, 1: 0},
        'col2': {0: 1, 1: 1},
        'col3': {0: 0, 1: 0},
        'total': {0: 6, 1: 2}
        })

@pytest.fixture
def fractions():
    return pd.DataFrame({
        'col1_match_fraction': {0: 1.0, 1: 0.0},
        'col2_match_fraction': {0: 1.0, 1: 1.0},
        'col3_match_fraction': {0: 0, 1: 0},
        'total_match_fraction': {0: 1.0, 1: 0.3333333333333333}
        })

@pytest.fixture
def cluster_nodes():
    cluster_nodes = pd.DataFrame({
        'id': {
            'col1_prot1': 'prot1', 'col1_prot3': 'prot3',
            'col2_prot1': 'prot1'},
        'fraction_clustered': {
            'col1_prot1': 1.0,
            'col1_prot3': 0.5,
            'col2_prot1': 1.0},
        'tag': {
            'col1_prot1': 'col1', 'col1_prot3': 'col1',
            'col2_prot1': 'col2'},
        'clust_id': {
            'col1_prot1': 0, 'col1_prot3': 0, 'col2_prot1': 0}
    })
    cluster_nodes.index.set_names('tagged_id',inplace=True)
    return cluster_nodes

def test_normalise_scores(score_list):
    expected = [
        pd.Series({
            'prot1':0.675,
            'prot2':0.600,
            'prot3':0.525,
        }),
        pd.Series({
            'prot4':0.75,
            'prot5':0.45,
            'prot6':0.60,
        })
    ]
    normed = mcl.normalise_scores(score_list)
    for i,normscores in enumerate(normed):
        assert_series_equal(normscores,expected[i])

def test_normalise_combined_scores():
    pass

def test_create_combined_network():
    pass

def test_run_MCL():
    pass

def test_separate_subclusters(clusters,nested_tags,separated_clusters):
    expected = separated_clusters
    res =  mcl.separate_subclusters(clusters,nested_tags)
    assert res == expected

def test_invert_nested_tags(nested_tags,inverted_tags):
    expected = inverted_tags
    res = mcl.invert_nested_tags(nested_tags)
    assert  res == expected 

def test_split_clusters(clusters,nested_tags,sample_split_clusters):
    expected = sample_split_clusters
    sample_tags = mcl.get_sample_tags(nested_tags)
    res = mcl.split_clusters(clusters,sample_tags)
    assert res == expected

def test_get_match_counts(nested_tags,mappings,sample_split_clusters,
                          counts,fractions):

    res_counts,res_fractions = mcl.get_match_counts(
                sample_split_clusters,mappings,nested_tags)

    assert_frame_equal(counts,res_counts)
    assert_frame_equal(fractions,res_fractions)

def test_get_cluster_nodes(separated_clusters,cluster_nodes):

    expected = cluster_nodes
    res = mcl.get_cluster_nodes(0,separated_clusters)

    assert_frame_equal(res,expected)
