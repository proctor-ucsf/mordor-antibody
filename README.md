# mordor-antibody
Analysis of multiplex antibody endpoints in the MORDOR Niger trial


## Description

This repository includes R code to run all of the analysis for the paper:

_Effect of biannual azithromycin distribution on antibody responses to malaria, bacterial, and protozoan pathogens in Niger_

currently in press at _Nature Communications_

A preprint of the article is available through medRxiv:

https://www.medrxiv.org/content/10.1101/2021.04.23.21255957v1

Should you have any questions about the files in this repository, please contact Ben Arnold at UCSF (ben.arnold@ucsf.edu).

## Linked Repositories and Additional Resources

### Open Science Framework
This GitHub repository is mirrored on the Open Science Framework (OSF).  The OSF project page includes additional study-related resources, including the pre-analysis plan, the datasets, and the compiled HTML computational notebooks created from the `.Rmd` files:

https://osf.io/954bt/

### Dryad 

The data have also been archived on Dryad: https://datadryad.org/stash/dataset/doi:10.7272/Q6VX0DSD

Arnold et al. (2021), Multiplex IgG antibody response and malaria parasitemia among children ages 1-59 months in the MORDOR Niger trial, 2015-2018, Dryad, Dataset, https://doi.org/10.7272/Q6VX0DSD

## _Nature_ Research Code Submission Items

Following: https://www.nature.com/documents/nr-software-policy.pdf

### System Requirements

All analyses were run using R software version 4.1.1 on Mac OSX Big Sur using the RStudio IDE (https://www.rstudio.com).

`> sessionInfo()`

`R version 4.1.1 (2021-08-10)`

`Platform: x86_64-apple-darwin17.0 (64-bit)`

`Running under: macOS Big Sur 11.6.1`

### Installation Guide

You can download and install R from CRAN: https://cran.r-project.org

You can download and install RStudio from their website: https://www.rstudio.com

All R packages required to run the analyses are sourced in the file `mordor-ab-Config.R`.

The installation time should be < 10 minutes total on a typical desktop computer.

### Instructions for Use

To reproduce all analyses in the paper, we recommend that you: 

1. clone the GitHub repository

2. Create a `data` subdirectory and copy the two datasets from OSF or Dryad

3. Create a `figures` subdirectory to store output. 

4. All of the analysis scripts should run smoothly (scripts `03-xx` to `16-xx`). 

The first two data processing scripts will not run because they read in our internal datasets with personally identifiable information, and create the final analysis datasets (shared publicly). We have included the data processing scripts for transparency and completeness.

You can run the `.Rmd` notebook scripts one-by-one or you can compile `mordor-ab-run-all.R`, which is the file we used to run the final analyses (e.g., from the command line `R CMD BATCH mordor-ab-run-all.R &`).

Running the all analyses on the above Mac desktop configuration required 31 minutes. 

Note that the only script that takes very long is `08-mordor-ab-community-means.Rmd` because estimating the ICCs and 95% CIs by bootstrapping binomial mixed models is computationally slow. 

Also note: we attempted to create a Binder virtual machine option for this project, but the underlying `.rds` dataset is so large (30 MB) that it took too long to spawn a remote docker container on the Binder server, and we didn't have time to troubleshoot/optimize.

### License

This project is covered by the CC0 1.0 Universal license.

