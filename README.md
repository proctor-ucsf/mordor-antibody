# mordor-antibody
Analysis of multiplex antibody endpoints in the MORDOR Niger trial


## Description

This repository includes R code to run all of the analysis for the paper:

_Effect of biannual azithromycin distribution on antibody responses to malaria, bacterial, and protozoan pathogens among children: A cluster-randomized, placebo-controlled trial in Niger_

(currently under peer review)

A preprint of the article is available through medRxiv:

https://www.medrxiv.org/content/10.1101/2021.04.23.21255957v1

This GitHub repository is mirrored on the Open Science Framework (OSF).  The OSF project page includes additional study-related resources, including the pre-analysis plan, the datasets, and the compiled HTML computational notebooks created from the `.Rmd` files:

https://osf.io/954bt/

To run the scripts: (1) clone the repo; (2) create a `data` subdirectory and copy the two datasets from the OSF project into it; (3) then, create a `figures` subdirectory to store output. All of the analysis scripts should run smoothly (scripts `03-xx` to `13-xx`). The first two data processing scripts will not run because they read in our internal datasets with PII. We have included them for transparency and completeness.

Should you have any questions about the files in this repository, please contact Ben Arnold at UCSF (ben.arnold@ucsf.edu).
