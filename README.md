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

CompaCt performs automated integrative comparative analysis of large-scale (protein) interaction datasets, identifying groups of interactors (e.g., protein complexes) in parallel in multiple species, allowing systematic identification of conserved as well as taxon-specific interactions.


...

For a more complete description of the software and its applications, please refer to the !manuscript!.

## Installation

CompaCt is implemented as a user-friendly command-line tool, as well as a flexible python package. Both are available after installing the python package and its dependencies. Installing with conda or using docker ensures these are available without requiring manual installation.

**dependencies**
- Python 3.6 or higher (check if that is correct, link..)
- rbo (link..)
- MCL (link..)

### pip
    pip install compact

### conda
    ???conda -c bioconda install compact

### docker
    docker pull joerivs/compact

# Usage

## Input Data

### interaction data

correlation, interaction, complexome profiling data

!!!! example interaction file


### orthology / mapping data

to enable comparison between species, or other datasets with different identifiers, pairwise orthology or identifier mappings need to be provided.

!!!! example mapping file


### settings file

To prevent an extensive list of input commands when analysing a large number of datasets, the list and path of all input files are provided in a settings file. This should be a plain text file, where each input file is listed per row, and fields are delimited by tab (\t) characters.

!!!!Description of format and fields

!!!!!example settings file


## Command line Tool

### when installed with pip or conda
    compact [options] path/to/settings.tsv

### when using docker
    ????docker run compact [options] path/to/settings.tsv

### command line arguments



## Python Package

    import compact ....
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

## Citing Compact
