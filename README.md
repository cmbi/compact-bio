# CompaCt

![Tests](https://github.com/joerivstrien/compact/actions/workflows/tests.yml/badge.svg)

### Comparative Clustering of (protein) interaction data

## Contents

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [License](./LICENSE)
- [Issues](#issues)
- [Citing CompaCt](#citing-compact)


## About

CompaCt performs automated integrative comparative analysis of large-scale (protein) interaction datasets, identifying groups of interactors (e.g., protein complexes) in parallel in multiple species, allowing systematic identification of conserved as well as taxon-specific interactions. For a more complete description of the software and its applications, please refer to the manuscript (see [here](#citing-compact)).

## Installation

CompaCt is implemented as a user-friendly command-line tool, as well as a flexible python package. Both are available after installing the python package and its dependencies. Installing with conda or using docker ensures the required dependencies are available without requiring manual installation.

**dependencies**
- [Python >= 3.7](https://www.python.org/) (check if that is correct, link..)
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
    git clone git@github.com:joerivstrien/compact.git
    cd compact
    pip install .

# Usage

## Input Data

### <ins> interaction data </ins>

CompaCt performs comparative clustering on a collection of (protein) correlation datasets. It expects a symmetric matrix of correlation or interaction scores, that represent interaction strength or likelihood. CompaCt does not require a specific metric or range of interaction scores, and can handle any numeric score as long as it can be used to sort or assign ranks to interactors.

The aim of CompaCt is to compare (protein) correlation datasets from different species or biological systems. To allow for prioritization of proteins that consistently cluster together in multiple experiments, multiple "replicated" correlation datasets that represent the same system can be provided. In the corresponding manuscript we have demonstrated that this can greatly improve quality of the CompaCt clustering results. A set of correlation datasets representing the same species or system ("replicates") are reffered to as a "collection". 

When using the command-line tool, the files containing these interaction scores should be provided as a tab-separated text file, with the first row and column containing the identifiers of the proteins,genes, etc. An example of an interaction file is available [here](LINK TO EXAMPLE interaction FILE).

Aside from directly providing correlation/interaction scores between proteins, CompaCt can automatically compute Pearson correlation scores from a matrix of expression/abundance data, like a complexome profile. These will then be used as interaction scores in subsequent CompaCt analysis. Similar to the interaction files, when using the command-line tool these data should be provided as a tab-separated text file, with the first row containing identifiers. An example of a file containing protein protein expression data is available [here](LINK TO EXAMPLE COMPLEXOME PROFILE)


### <ins>pairwise orthology data</ins>

to enable comparison between species, or other datasets with different identifiers, pairwise orthology or identifier mappings need to be provided. The identifier pairs should correspond to the identifiers used in the interaction or expression data. When using the command-line tool, the identifier pairs should be provided as a text file, with two columns separated by a tab. each column contains identifiers of one of the two species this orthology file relates to. An example mapping file is available [here](LINK TO EXAMPLE MAPPING FILE)

### <ins>settings file</ins>

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

The ordering of the files in the input settings file does not matter, as long as the fields are correctly describing it. An example settings file is available [here](LINK TO EXAMPLE SETTINGS FILE)

### <ins>annotation reference file</ins>

!!!!STILL TO DESCRIBE

## Command line Tool

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

### Output

!!! STILL TO ADD description of the output files etc. 


## Python Package

    import compact ....
    MAKE AND EXAMPLE SCRIPT DOING SOME BASIC STUFF. WITH COMMENTS
    .....

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
If you have questions or encounter any problems or bugs, please report this in the [issue channel](https://github.com/joerivstrien/compact/issues).


## Citing Compact

Publication Pending