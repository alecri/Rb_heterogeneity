### R script to reproduce the results in:
###   Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. 
###   A new measure of between-studies heterogeneity in meta-analysis. 
###   Statistics in Medicine. 2016 May 10. doi: 10.1002/sim.6980. 

## Packages required
biasl <- list()
coveragel <- list()
meanMeasurel <- list()
for (i in 1:3){
  load(paste0("data/170728_main_sim", i, ".RData"))
  biasl[[i]] <- biasSum
  coveragel[[i]] <- coverageSum
  meanMeasurel[[i]] <- meanMeasure
  rm(list = setdiff(ls(), c("biasl", "coveragel", "meanMeasurel")))
}

library(tidyverse)
x = coveragel[[1]]

lapply(biasl, function(x)
  x %>%
    select(-RR, -tau2) %>%
    rename(Rb = R, K = S) %>%
    select(Rb, K, everything()) %>%
    filter(Rb %in% c(.1, .5, .7), K %in% c(5, 20, 50, 100)) %>%
    round())

lapply(coveragel, function(x)
  x %>%
    select(-RR, -tau2) %>%
    rename(Rb = R, K = S) %>%
    select(Rb, K, everything()) %>%
    filter(Rb %in% c(.1, .5, .7), K %in% c(5, 20, 50, 100)) %>%
    round())

