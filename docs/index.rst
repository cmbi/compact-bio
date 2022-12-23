.. compact documentation master file, created by
   sphinx-quickstart on Fri Mar  4 10:32:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CompaCt's documentation!
===================================

`CompaCt`_ performs automated integrative comparative analysis of large-scale 
(protein) interaction datasets, identifying groups of interactors (e.g., protein complexes)
in parallel in multiple species, allowing systematic identification and comparison of 
conserved as well as taxon-specific components of protein complexes and other interactions.

Package modules
---------------

   - ``main``: high-level workflow of this project, to run the major steps of the CompaCt analsysis
   - ``pairwise_scoring``: perform pairwise rbo scoring between interaction datasets
   - ``reciprocal_hits``: find protein pairs that are reciprocal “top” hits after rbo scoring
   - ``MCL_clustering``: create combined network from rbo scored data, cluster with MCL, and processes cluster results
   - ``process_data``: general functions for processing,parsing used file formats and data sctructures.
   - ``utils``: module with generally useful functions used in multiple modules

.. _`CompaCt`: https://github.com/joerivstrien/compact


API Reference
-------------
.. toctree::
   :maxdepth: 2
   
   compact.main
   compact.pairwise_scoring
   compact.reciprocal_hits
   compact.MCL_clustering
   compact.process_data
   compact.utils

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
