# CompaCt

![Tests](https://github.com/joerivstrien/compact-bio/actions/workflows/tests.yml/badge.svg)

### Comparative Clustering of (protein) interaction data

## Contents

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Python Package](#python-package)
- [License](#licence)
- [Issues](#issues)
- [Citing CompaCt](#citing-compact)


## About

CompaCt performs automated integrative comparative analysis of large-scale (protein) interaction datasets, identifying groups of interactors (e.g., protein complexes) in parallel in multiple species, allowing systematic identification and comparison of conserved as well as taxon-specific components of protein complexes and other interactions. For a more complete description of the software and its applications, please refer to the manuscript (see [here](#citing-compact)).

## Installation

CompaCt is implemented as a user-friendly command-line tool, as well as a flexible python package. Both are available after installing the python package and its dependencies. Installing with conda or using docker ensures the required dependencies are available without requiring manual installation.

**dependencies**
- [Python >= 3.7](https://www.python.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [rbo](https://github.com/changyaochen/rbo)
- [MCL](https://micans.org/mcl/)
- [optional]: [fastrbo](https://github.com/joerivstrien/fastrbo) -- optional alternative to "rbo" package
    -  faster implementation of rank biased overlap, manual installation required. Follow the link for installation instructions. Fastrbo is automatically available when using the docker image

### pip
    publication on pip is pending

### conda
    publication on conda is pending

### docker
    docker image is pending

### installation from repository
    git clone git@github.com:joerivstrien/compact-bio.git
    cd compact-bio
    pip install .

# Usage

## Input Data

Below follows a description of the different types of input files that the CompaCt tool expects as input. The format of all input files is the tab-separated plain text format, as these are easy to create/edit with both spreadsheet software (e.g., excel, google-docs etc.) and programming languages like python or R.


### interaction data

CompaCt performs comparative clustering on a collection of (protein) correlation datasets. It expects a symmetric matrix of correlation or interaction scores, that represent interaction strength or likelihood. CompaCt does not require a specific metric or range of interaction scores, and can handle any numeric score as long as it can be used to sort or assign ranks to interactors.

The aim of CompaCt is to compare (protein) correlation datasets from different species or biological systems. To allow for prioritization of proteins that consistently cluster together in multiple experiments, multiple "replicated" correlation datasets that represent the same system can be provided. In the corresponding manuscript we have demonstrated that this can greatly improve quality of the CompaCt clustering results. A set of correlation datasets representing the same species or system ("replicates") are reffered to as a "collection". 

When using the command-line tool, the files containing these interaction scores should be provided as a tab-separated text file, with the first row and column containing the identifiers of the proteins,genes, etc. An example of an interaction file is available [here](https://github.com/joerivstrien/compact-bio/blob/master/example_data/correlated_profile.tsv).

Aside from directly providing correlation/interaction scores between proteins, CompaCt can automatically compute Pearson correlation scores from a matrix of expression/abundance data, like a complexome profile. These will then be used as interaction scores in subsequent CompaCt analysis. Similar to the interaction files, when using the command-line tool these data should be provided as a tab-separated text file, with the first row containing identifiers. An example of a file containing protein protein expression data is available [here](https://github.com/joerivstrien/compact-bio/blob/master/example_data/abundance_profile.tsv)


### pairwise orthology data

to enable comparison between species, or other datasets with different identifiers, pairwise orthology or identifier mappings need to be provided. The identifier pairs should correspond to the identifiers used in the interaction or expression data. When using the command-line tool, the identifier pairs should be provided as a text file, with two columns separated by a tab. each column contains identifiers of one of the two species this orthology file relates to. An example mapping file is available [here](https://github.com/joerivstrien/compact-bio/blob/master/example_data/pairwise_orthology.tsv)

### settings file

To prevent an extensive list of input commands when performing CompaCt analysis using the command-line tool on a large number of datasets, the list and path of all input files are provided in a settings file. This should be a "tab" delimited plain text file, where each input file is listed per row. Each file is desribed with four fields (delimited by tab characters). below follows a description of these four fields for interaction and orthology files.

<ins>interaction data file</ins>

1. file type identifier: should be "INT" specifying that this is an interaction data file.
2. "collection" identifier: an identifier for the "collection" (e.g., species, cell-type, tissue, etc.) this correlation dataset represents.
3. "replicate" identifier: an identifier for this specific interaction dataset.
4. the location and name of the file (e.g., path/to/file.tsv)   

<ins>pairwise orthology file</ins>

1. file type identifier: should be "ORTH" specifying that this is a pairwise orthology file.
2. FROM collection identifier: the collection identifier for the collection whose identifiers are in the first column of the pairwise orthology file.
3. TO collection identifier: the collection identifier for the collection whose identifiers are in the second column of the pairwise orthology file.
4. the location and name of the file (e.g., path/to/file.tsv)   

The ordering of the files in the input settings file does not matter, as long as the fields are correctly describing it. An example settings file is available [here](https://github.com/joerivstrien/compact-bio/blob/master/example_data/input_settings.tsv)

### <ins>annotation reference file</ins>

CompaCt is optionally able to automatically annotate the resulting clusters based on overlap with a provided reference set of protein complexes or pathways. The reference should be provided in a file with the [GMT format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29), a commonly used and available format to store gene sets for use in for example gene set enrichment analysis. The identifiers used in the reference file should match those of one of the analysed collections/species. An example GMT file with reference complexes is available [here](https://github.com/joerivstrien/compact-bio/blob/master/example_data/reference_complexes.gmt).

## CompaCt Command line Tool

### when installed with pip or conda
    compact [options] path/to/settings.tsv

### when using Docker
    coming soon

### Command Line Arguments
    #################################
    #### commonly useful options ####
    #################################

    -h, --help          show a complete list of all command line arguments
    -v, --version       show program's version
    -i {abun,int}, --in-type {abun,int}
                        whether provided input interaction data are abundance/expression
                        values or correlation matrices. default='abun'
    -o [OUTPUT_LOC], --output-loc [OUTPUT_LOC]
                        preferred location of the result folder. defaults to current
                        working directory
    -j JOB_NAME, --job-name JOB_NAME
                        name of current job, used in name of result folder name.
                        defaults to concatenated collection names
    -t [PROCESSES], --processes [PROCESSES]
                        number of threads or processes to use. default=1
    -p [P]              rank biased overlap p parameter, value between 0 and 1. 
                        default=0.9.
                        determines "top-heaviness" of the rank biased overlap metric.
                        Lower values correspond to more weight towards the top of the
                        ranked interactor profilewhen computing similarity.
    --mcl-inflation [MCL_INFLATION]
                        inflation param of MCL tool, determines cluster granularity. default=2. 
                        higher values will typically result in more coarse-grained clusters
    --perf-cluster-annot
                        perform automatic annotation of clusters using reference
                        groups
    --ref-file [REF_FILE]  
                        gmt format file with reference groups, when perf-cluster-
                        annot is used
    --ref-tag [REF_TAG]   tag of collection that the annotation reference represents.
    

    ##############################################################
    #### more obscure options, ignoring these is usually fine ####
    ##############################################################

    -m [MIN_SEARCH_WEIGHT], --min-search-weight [MIN_SEARCH_WEIGHT]
                        minimum search weight for rank biased overlap, default=0.99
    --th-criterium {percent,best}
                        threshold criterium for determination of reciprocal top hits.
                        default='percent'
    --th-percent [TH_PERCENT]
                        if th-criterium='percent': top % to consider top hits.
                        default=1.
    --include_within    include within-sample interaction scores in clustering input.
                        default=False
    --wbratio [WBRATIO]   ratio of within/between score averages. default=1
    --skip-filter-clusters
                        use if low-scoring clusters should not be filtered
    --annot-ref-th [ANNOT_REF_TH]
                        min fraction of reference group that must be present in
                        cluster. default=0.5
    --annot-mem-th [ANNOT_MEM_TH]
                        min fraction-present of clust members to be counted for
                        annot-ref-th. default=0.25

## Output

When running the compact command line tool, the results will be saved to a results folder containing several files. We will describe the content of each file in more detail below.
All result files are in tab separated text format (.tsv), so they are easy to load into a spreadsheet software (e.g., google docs, excel, etc.) or work with using a programming language like R or python. 

#### <ins> clust info file</ins>: "clust_info.tsv"

This file contains the full list of result clusters that passed CompaCt's filtering steps.
It contains various columns with information regarding the cluster, like their size, coherence and other metrics. Below follows a description of the information available in each of the various columns:
- *_size: number of members of each subcluster and total cluster
- n_represented: number of complexomes that have cluster members in this cluster
- *_over_0.5_size: number of members with a fraction clustered score (FC) of 0.5 or greater
- robust_represented: number of complexomes that have cluster members in this cluster with a FC of 0.5 or greater
- *_match_fraction: the fraction of total possible “matches” within each complexome, or overall. A metric to quantify the coherence of a cluster across datasets

#### <ins>cluster member files</ins>: "*_cluster_members.tsv"
For each collection of datasets (e.g., all datasets corresponding to the same species) that was provided as input, a cluster member file is available, containing all the clustered proteins and their cluster assignments, along with additional information like their orthologs in other species as well as their FC (fraction clustered) scores. Below follows a description of the information available in each of the various columns:

- clust_id: id of cluster this protein/gene entry is clustered with
- id: identifier of clustered protein/gene
- fraction_clustered: the FC score of this protein with this cluster. Defined as: number of replicate datasets in which this protein clusters with this cluster divided by the total number of replicates for this collection
- *_mapping: orthologous protein in other collection, based on provided pairwise orthology
- best_guess_selection: whether this protein is part of the "best guess" selection criterion of the respective cluster
- match_over_threshold: whether this protein has an ortholog in another collection that has a FC score over 0.5

#### <ins>MCL result file</ins>: "MCL_result.tsv"
The raw cluster output from the MCL clustering software. for each line contains the members of one of the MCL clusters, separated by tab characters. The protein/gene identifiers used are prepended by the replicate dataset identifier that they originate from. for more details regarding the MCL clustering software and its output please refer to [their documentation](https://micans.org/mcl/)

#### <ins>combined network file</ins>: "combined_network.tsv"
The combined network generated by compact, containing proteins/genes from all provided datasets, to be used as input for clustering by MCL. The protein/gene identifiers used are prepended by the replicate dataset identifier that they originate from.

#### <ins>clust_nodes</ins>: "clust_nodes.tsv"
contains all clustered proteins, along with additional information. Can be used together with the clust_eges.tsv file for visualization with tools like cytoscape or analysis of CompaCt's output clusters

#### <ins>clust_edges</ins>: "clust_edges.tsv"
contains all edges between proteins part of the same cluster, along with additional information. Can be used together with the clust_nodes.tsv file for network vizualisation with tools like cytoscape or analysis of CompaCt's output clusters 

## Python Package

Aside from using the CompaCt command line tool, that performs the complete analysis from expression/abundance datasets to annotated species-specific clusters, CompaCt can also be used as a python package. From python CompaCt can be used more flexibly, for example to rerun only specific steps in the analysis with changed parameters. For complete documentation of all the modules and functions available in the CompaCt package refer to the [package documentation](LINK TO PACKAGE DOCUMENTATION)

Running a complete CompaCt analysis from python

    from compact.main import main
    import compact.process_data as prd
    import compact.utils as ut

    # parse complexome profiles and orthology files as listed in settings.tsv
    sample_data,mapping_data = prd.parse_settings('path/to/settings_file.tsv')
    samples = prd.parse_profiles(sample_data)
    mappings = prd.parse_mappings(mapping_data)

    # get collection/replicate structure from settings file
    nested_tags = prd.get_nested_tags(sample_data)

    # generate Pearson correlation matrices from the complexome profile abundance data
    corr_matrices = ut.correlate_samples(samples)

    # run complete CompaCt analysis
    main(
        nested_tags,corr_matrices,mappings,
        output_location='results/',
        job_name='example_run',
        processes=4)

Now let's say you would like to use a different granularity for the clustering step of the analysis, without having to rerun the computationally expensive rbo score computation. You could do this as follows:

    from compact.MCL_clustering import run_MCL
    from compact.main import process_mcl_result,save_results

    # rerun MCL with a different inflation parameter value
    run_MCL('results/example_run_results/combined_network.tsv',
            'results/rerun_results/mcl_rerun_result.tsv',
             inflation=2.5,processes=4)

    # process the MCL results
    # to get scored and filtered species-specific clusters
    mcl_res = process_mcl_result(
        'results/rerun_results/mcl_rerun_result.tsv',nested_tags,
        'results/example_run_results/combined_network.tsv',mappings)

    # save the results to an output folder for manual inspection
    save_results(mcl_res,'../results/rerun_results/',mappings)

## Licence

    CompaCt -- Comparative Clustering
    Copyright (C) 2022 Radboud University Medical Center

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Issues
If you have questions or encounter any problems or bugs, please report them in the [issue channel](https://github.com/joerivstrien/compact-bio/issues).


## Citing Compact

Publication Pending