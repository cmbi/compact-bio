from compact import pairwise_scoring as ps
import pytest

@pytest.fixture
def data():
    return {
        'left':['a','b','c'],
        'right': ['b','c','d'],
        'mapping':{'a':'b'},
    }



def test_rename_indices(data):
    expected = (
        ['a|b','b','c'],
        ['a|b','c','d']
        )
    assert ps.rename_indices(
        data['left'],data['right'],data['mapping']) == expected
    
def assign_matched():
    pass

def score_comparison():
    pass
