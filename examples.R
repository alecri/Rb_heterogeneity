### R script to reproduce the results in:
###   Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. 
###   A new measure of between-studies heterogeneity in meta-analysis. 
###   Statistics in Medicine. 2016 May 10. doi: 10.1002/sim.6980. 

## Packages required
#devtools::install_github("alecri/hetmeta")
library(metafor)
library(hetmeta)


## -----------------------------------------------------------------------------
## Example 5.1
## Self-management education and regular medical review for adults with asthma

data(dat.gibson2002)
#help(dat.gibson2002)

dat_gibson2002 <- escalc(measure = "SMD", m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, 
              sd2i = sd2i, n2i = n2i, data = dat.gibson2002)
res_gibson2002 <- rma(yi = yi, vi = vi, data = dat_gibson2002, method = "DL")
res_gibson2002

het_gibson2002 <- hetmeta(res_gibson2002)
het_gibson2002
confint(het_gibson2002, rma.type = TRUE, level = 95)


## -----------------------------------------------------------------------------
## Example 5.2
## Bacillus Calmette-Guerin vaccine and turberculosis

data(dat.bcg)
#help(dat.bcg)

dat_bcg <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, 
                  data = dat.bcg)
res_bcg <- rma(yi, vi, data = dat_bcg, method = "DL")
res_bcg
predict(res_bcg, transf = exp)

het_bcg <- hetmeta(res_bcg)
het_bcg
confint(het_bcg, rma.type = T, level = 95)



## -----------------------------------------------------------------------------
## Example 5.3
## Circumcision status and risk of sexually transmitted infections

load("data/millett.Rdata")
millett

res_millett <- rma(yi = yi, vi = vi, data = millett, method = "DL")
res_millett
predict(res_millett, transf = exp)


het_millett <- hetmeta(res_millett)
het_millett
confint(het_millett, rma.type = F, level = 95)


## -----------------------------------------------------------------------------
## Table I
## Example of a hypothetical meta-analysis of 10 studies.

## function for typical s2 as defined in I2
s2typ <- function(si2){
  k <- length(si2)
  wi <- 1/si2
  (k-1)*sum(wi)/(sum(wi)^2 - sum(wi^2))   
}

A <- c(6, 6.1, 6.2, 5.9, 6.0, 5.9, 6.1, 5.8, 6.0, 6.2)
B <- c(5, 19, 3, 15, 6, 23, 4, 17, 2, 8.8)

round(data.frame(
  CV_vi = c(sd(A)/mean(A), sd(B)/mean(B)),
  s2_1 = c(s2typ(A), s2typ(B)),
  s2_2 = c(length(A)/(sum(1/A)), length(B)/(sum(1/B))),
  v_mean = c(mean(A), mean(B))
), 3)
