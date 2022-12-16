from compact import pairwise_scoring as ps


def test_rename_indices():
    left = ['a','b','c']
    right = ['b','c','d']
    mapping = {'a':'b'}
    expected = (
        ['a|b','b','c'],
        ['a|b','c','d']
        )
    assert ps.rename_indices(
        left,right,mapping) == expected