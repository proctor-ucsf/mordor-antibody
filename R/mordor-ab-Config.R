
#----------------------------
# mordor-ab-Config.R
#
# Analysis of multiplex antibody 
# endpoints in the MORDOR Niger trial
#
# Configuration file
#----------------------------


#----------------------------
# load packages
#----------------------------

library(here)
library(kableExtra)
library(tidyverse)
library(gridExtra)
library(viridisLite)

#----------------------------
# stats packages
#----------------------------
library(mixtools)
library(mgcv)
library(coin)

#----------------------------
# set up for parallel computing
#----------------------------
library(foreach)
library(doParallel)
registerDoParallel(detectCores() - 1)


#----------------------------
# custom color pallettes
#----------------------------

# safe color blind palette
# http://jfly.iam.u-tokyo.ac.jp/color/
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# Reference: Bang Wong, Nature Methods 2011: https://www.nature.com/articles/nmeth.1618
cbpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Viridis color palette
vircols <- viridis(9,begin=0,end=0.97)
