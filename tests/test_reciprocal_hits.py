from compact import reciprocal_hits as rh
import pytest
import pandas as pd
from pandas.testing import assert_series_equal

@pytest.fixture
def scores():
    return pd.DataFrame({
        'protB1':{'protA1':0.9,'protA2':0.3,'protA3':0.5},
        'protB2':{'protA1':0.6,'protA2':0.9,'protA3':0},
        'protB3':{'protA1':0.5,'protA2':0,'protA3':0.9}
    })

@pytest.fixture
def top_hits():
    return pd.Series({
        ('protA1', 'protB1'): 0.9,
        ('protA2', 'protB2'): 0.9,
        ('protA3', 'protB3'): 0.9        
    })    

def test_get_reciprocal_top_hits(scores,top_hits):
    expected = top_hits
    res = rh.get_reciprocal_top_hits(scores,percent=30)
    assert_series_equal(res,expected)


