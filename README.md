## Overview

This repository reproduces the R analyses from [Tsakiroglou et al. 2024](https://doi.org/10.3390/ijms252312653).

## Installation

Packages used in the analyses can be installed by running `install/install_r_libs.R` in a fresh R session.

The package versions used for the manuscript analyses can be found within `output/sessionInfo` for each script.

## Directions

Scripts to run the analyses from the manuscript can be found in the `code` folder. Outputs from each analysis are already included in this repository, so each script can be run as a standalone analysis.

To reproduce the analysis from scratch using the data files in `raw`, each script should be run in the order they are numbered. Raw CEL files can be obtained manually from Gene Expression Omnibus (GEO) accession no. GSE276395 and saved in `raw/GSE276395_RAW`, or automatically using the R script `code/01_process_cel.R`.

[![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg