
#-----------------------------------
# This script runs all computational
# notebooks used in the article:
#
# Arzika et al.
# Biannual azithromycin distribution 
# and antibody responses to malaria, bacterial, 
# and protozoan pathogens among Nigerien children
#
# There are 16 scripts, run sequentially
# and sequentially numbered in the repository
#
# NOTE: The first 2 scripts cannot be
# run on the publicly available data
# They create the public datasets, which
# are used for all analyses. 
# We have included them and their
# output files (.html) for completeness.
#-----------------------------------
library(here)
here::here()

#-----------------------------------
# Raw study data processing
# to create public data files
#
# (not run for public replication)
#-----------------------------------
rmarkdown::render(here::here("R/01-mordor-ab-format-data-gps-climate.Rmd"),
                  output_file = here::here("R/01-mordor-ab-format-data-gps-climate.html"))


rmarkdown::render(here::here("R/02-mordor-ab-format-data.Rmd"),
                  output_file = here::here("R/02-mordor-ab-format-data.html"))

#-----------------------------------
# Figure 1. 
# CONSORT flow
#-----------------------------------
rmarkdown::render(here::here("R/03-mordor-ab-consort.Rmd"),
                  output_file = here::here("R/03-mordor-ab-consort.html"))

#-----------------------------------
# Table 1. 
# Baseline balance
#-----------------------------------
rmarkdown::render(here::here("R/04-mordor-ab-baseline-balance.Rmd"),
                  output_file = here::here("R/04-mordor-ab-baseline-balance.html"))

#-----------------------------------
# SI Figures 1 and 2
# descriptive summaries
#-----------------------------------
rmarkdown::render(here::here("R/05-mordor-ab-descriptive-summary.Rmd"),
                  output_file = here::here("R/05-mordor-ab-descriptive-summary.html"))

#-----------------------------------
# SI Figure 3
# measurement timing
#-----------------------------------
rmarkdown::render(here::here("R/06-mordor-ab-measurement-timing.Rmd"),
                  output_file = here::here("R/06-mordor-ab-measurement-timing.html"))

#-----------------------------------
# SI Figure 4
# pairwise correlation of IgG responses
#-----------------------------------
rmarkdown::render(here::here("R/07-mordor-ab-pairwise-ab-correlation.Rmd"),
                  output_file = here::here("R/07-mordor-ab-pairwise-ab-correlation.html"))

#-----------------------------------
# Community level means
#-----------------------------------
rmarkdown::render(here::here("R/08-mordor-ab-community-means.Rmd"),
                  output_file = here::here("R/08-mordor-ab-community-means.html"))

#-----------------------------------
# Age dependent antibody curves
#-----------------------------------
rmarkdown::render(here::here("R/09-mordor-ab-age-curves.Rmd"),
                  output_file = here::here("R/09-mordor-ab-age-curves.html"))

#-----------------------------------
# Primary analysis
#-----------------------------------
rmarkdown::render(here::here("R/10-mordor-ab-mean-diffs.Rmd"),
                  output_file = here::here("R/10-mordor-ab-mean-diffs.html"))

rmarkdown::render(here::here("R/11-mordor-ab-malaria-composite.Rmd"),
                  output_file = here::here("R/11-mordor-ab-malaria-composite.html"))

#-----------------------------------
# Longitudinal analyses
#-----------------------------------
rmarkdown::render(here::here("R/12-mordor-ab-longitudinal.Rmd"),
                  output_file = here::here("R/12-mordor-ab-longitudinal.html"))

#-----------------------------------
# Pre-specified Subgroup analyses
#-----------------------------------
rmarkdown::render(here::here("R/13-mordor-ab-subgroup-analyses.Rmd"),
                  output_file = here::here("R/13-mordor-ab-subgroup-analyses.html"))

#-----------------------------------
# Sensitivity analysis that varies
# MFI cutoff values by +/- 20%
#-----------------------------------
rmarkdown::render(here::here("R/14-mordor-ab-cutoff-sensitivity-analysis.Rmd"),
                  output_file = here::here("R/14-mordor-ab-cutoff-sensitivity-analysis.html"))

#-----------------------------------
# Sensitivity analysis that conducts
# a leave-one-out (LOO) analysis
# and computes jackknife estimates
# of the bias
#-----------------------------------
rmarkdown::render(here::here("R/15-mordor-ab-loo-sensitivity-analysis.Rmd"),
                  output_file = here::here("R/15-mordor-ab-loo-sensitivity-analysis.html"))

#-----------------------------------
# Comparison of malaria parasitemia
# and P. falciparum seroprevalence
# and (not-prespecified) subgroup
# analysis by baseline malaria
# parasitemia
#-----------------------------------
rmarkdown::render(here::here("R/16-mordor-ab-pf-infection-comparison.Rmd"),
                  output_file = here::here("R/16-mordor-ab-pf-infection-comparison.html"))



