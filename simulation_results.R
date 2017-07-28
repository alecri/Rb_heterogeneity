### R script to reproduce the results in:
###   Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. 
###   A new measure of between-studies heterogeneity in meta-analysis. 
###   Statistics in Medicine. 2016 May 10. doi: 10.1002/sim.6980. 

## Packages required
library(tidyverse)


## Loading results of the simulation study design
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

## Table II (percent relative bias)
lapply(biasl, function(x)
  x %>%
    select(-RR, -tau2) %>%
    rename(Rb = R, K = S) %>%
    select(Rb, K, everything()) %>%
    filter(Rb %in% c(.1, .5, .7), K %in% c(5, 20, 50, 100)) %>%
    round())

## Table III (confidence interval coverage)
lapply(coveragel, function(x)
  x %>%
    select(-RR, -tau2) %>%
    rename(Rb = R, K = S) %>%
    select(Rb, K, everything()) %>%
    filter(Rb %in% c(.1, .5, .7), K %in% c(5, 20, 50, 100)) %>%
    round())


## html tables
library(htmlTable)
tableListBias <- lapply(biasl, function(x) 
  htmlTable(round(x[, -c(1, 3, 4)], 1),
            cgroup = c("", "CV = .5", "CV = 1", "CV = 2"), 
            n.cgroup = c(1, 3, 3, 3), align = 'rccccccccc',
            header = c("S", rep(c("Rb", "I2", "Ri"), 3)),
            rgroup = c("R = .1", "R = .3", "R = .5", "R = .7", "R = .9"),
            n.rgroup = c(rep(6, 5)),
            rnames = rep("", nrow(x))       
  )
)
interactive()
lapply(tableListBias, print)