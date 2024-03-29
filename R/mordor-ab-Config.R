
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
library(lubridate)
library(table1)
library(gridExtra)
library(viridisLite)
library(RColorBrewer)

#----------------------------
# stats packages
#----------------------------
library(mixtools)
library(mgcv)
library(coin)
library(rptR)

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

# NY Times rainbow palettes from the Upshot
# https://www.nytimes.com/interactive/2020/03/21/upshot/coronavirus-deaths-by-country.html
nytpal <- c("#510000", "#AC112D", "#EC6D47", "#F2A058", "#F7D269", "#839772", "#325D8A")

# https://www.nytimes.com/2020/04/23/us/coronavirus-early-outbreaks-cities.html
nytpal2 <- c("#7486AD", "#77658E","#A1493B","#C17B3F", "#E6CA83")

# Viridis color palette
vircols <- viridis(9,begin=0,end=0.97)
