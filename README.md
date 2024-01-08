# Clonalscope

Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description

Clonalscope is a subclone detection method based on copy number alterations (CNAs) for single-cell and ST tumor sequencing data. Clonalscope leverages prior information from matched bulk DNA-seq data, and implements a nested Chinese Restaurant Process to model the evolutionary process in tumors. Clonalscope is able to detect subclones, label malignant cells (with the help of matched WGS/WES data), and trace subclones for both scRNA-seq and ST data.

For more information about the method, please check out the [manuscript](https://doi.org/10.1101/2022.07.05.498882). <br/>

## Overview of Clonalscope

![](inst/plots/overview.png?raw=true "Overview of subclone detection with Clonalscope")

## Install

-   You can install Clonalscope with the code below:

``` r
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install.packages("devtools")
options(timeout=9999999)
devtools::install_github("seasoncloud/Clonalscope") # install
library(Clonalscope) # load

# It takes <5 mins for the installation.
```

-   You can download example datasets for Clonalscope with the following command:

    Using terminal, download the repository.

        git clone https://github.com/seasoncloud/Clonalscope.git

-   We have included example data in the folder data-raw/.

## Detailed tutorials with example datasets

-   Click the links below to read detailed tutorials for different data types.

1.  [scRNA-seq subclone detection and malignant cell labeling](https://github.com/seasoncloud/Clonalscope/tree/main/samples/P5931/scRNA)
2.  [Spatial Transcriptomic subclone detection](https://github.com/seasoncloud/Clonalscope/tree/main/samples/BC_ductal2/ST)
3.  [Spatial Transcriptomic subclone tracing](https://github.com/seasoncloud/Clonalscope/tree/main/samples/BC_ductal1/ST)
4.  [scRNA-seq subclone detection and malignant cell labeling without paired WGS/WES](https://github.com/seasoncloud/Clonalscope/tree/no_pairedWGS_branch/samples/P5847/scRNA) <br/>

## Citation

Wu, C.-Y. et al. Cancer subclone detection based on DNA copy number in single cell and spatial omic sequencing data. *bioRxiv* (2022): <https://doi.org/10.1101/2022.07.05.498882>
