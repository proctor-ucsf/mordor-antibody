##################################################
# MORDOR
# Create DBS Dataset
# Victoria Le
# Dec 2019
##################################################

setwd("~/Box Sync/MORDOR/MORDOR Morbidity Serology Dataset")

library(tidyverse)

morbidity <- read_csv("data/MORDOR-Morbidity-merged-20Mar2019.csv",
                      col_types = cols(arm = col_character(),
                                       txlabelgwu = col_character()))

# Select relevant variables and drop obs without bspot code
vars <- c("masterphase", "masterperson_mcd", "code_bspot_", "gwu", "gwunotes_", "arm", "txlabelgwu", 
          "studypop_", "phaseidmep", "dobmcdagemos", "agecc_mcd", "gendermcd_", "thicksmear_cons_bsmear_",
          "officialdensity_bsmear_", "gameto_bsmear_")

morbiditydbs <- morbidity %>% 
                select(vars) %>% 
                filter(is.na(gwunotes_)) %>% 
                mutate(gwu = str_to_upper(gwu))

# Load census data set (includes vital status and time at risk)
census <- read_csv("data/Stata-mordor-Niger-long-dedupe-forYing.csv")

# Combine masterperson and phase to get masterphase
census1 <- census %>% 
           mutate(masterphase = str_c(census$masterperson, census$phase, sep = ".")) %>% 
           select(masterphase, statusx, timeatriskx)

# Merge morbidity and census data by masterphase
dbs <- left_join(morbiditydbs, census1, by = "masterphase")

dbs <- arrange(dbs, masterperson_mcd, phaseidmep)

write_csv(dbs, "MORDOR-morbidity-serology-2019-12-06.csv")
