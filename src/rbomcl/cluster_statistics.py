"""
module to perform statistics on clusters

to be used on raw clusters, before aggregating replicates etc.


given a cluster of a certain composition
    - determine number of matches within this cluster
        - code to count total matches
        - requires mappings, nested tags
        - take inspiration from process_cluster code that does something similar
    - sample a number of clusters randomly with the same distribution over samples
    - determine null distribution of number of matches
    - compute z-score, p-value
        - what fraction of nulls >= the scores? (p-value)
            - this doesn't make any assumptions about the distribution of n_matches
    - Then you need to do multiple testing correction to get FDR
        - the code to convert p-vals into q-vals is already complete in ../Apicomplexa_project/code/significance.py. Take from there!
            - CREATE SEPARATE MODULE FOR (MP) FDR CORRECTION. CAN BE USEFUL IN MANY CASES!!!!


SO. STEPS:
    - Determine what exact input you need, in what structure.
    - Then go from there
"""




