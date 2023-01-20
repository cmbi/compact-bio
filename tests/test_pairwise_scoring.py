from compact import pairwise_scoring as ps
import pytest

@pytest.fixture
def data():
    return {
        'left':['a','b','c'],
        'right': ['b','c','d'],
        'mapping':{'a':'b'},
    }

@pytest.fixture
def comparison():
    return (
        ('prot1',['A','B','C']),
        ('prot2',['B','A','D'])
    )

def test_rename_indices(data):
    expected = (
        ['a|b','b','c'],
        ['a|b','c','d']
        )
    assert ps.rename_indices(
        data['left'],data['right'],data['mapping']) == expected
   

def test_score_comparison(comparison):
    expected = (('prot1','prot2'),0.14399999999999996)
    assert ps.score_comparison(comparison,0.9) == expected
