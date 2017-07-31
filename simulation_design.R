### R script to reproduce the results in:
###   Crippa A, Khudyakov P, Wang M, Orsini N, Spiegelman D. 
###   A new measure of between-studies heterogeneity in meta-analysis. 
###   Statistics in Medicine. 2016 May 10. doi: 10.1002/sim.6980. 

## Packages required
library(reshape)


## -----------------------------------------------------------------------------
## Functions used in the simulation/analysis

## Alternative function to conduct meta-analysis
het_measures <- function(yi, vi){
  si <- function(n) sum(wfi^n)
  k <- length(yi)
  wfi <- 1/vi
  y_fixed <- sum(yi*wfi)/sum(wfi)
  ##Q-test and DL tau2 estimates
  Q <- sum(wfi*(yi-y_fixed)^2)
  tau2 <- max((Q - (k-1))/(sum(wfi)- sum(wfi^2)/sum(wfi)), 0)
  wri <- 1/(vi + tau2)
  y_random <- sum(yi*wri)/sum(wri)
  se_yr <- sqrt(1/sum(wri))
  
  ## Measures of heterogeneity
  Rb <- tau2/(k*se_yr^2)
  I2 <- max((Q-(k-1))/Q, 0)
  Ri <- tau2/(tau2 + k/si(1))
  CVb <- sqrt(tau2/se_yr^2)
  cv2 <- k*si(2)/(si(1)^2) - 1
  s2I2 <- (k-1)*si(1)/(si(1)^2 - si(2))
  s2Ri <- k/si(1)
  
  ## Approximate asymptotic stderr
  se.Q <- sqrt(2*(k-1) + 4*(si(1) - si(2)/si(1))*tau2 + 2*(si(2) - 2*si(3)/si(1) + si(2)^2/si(1)^2)*tau2^2)
  se.tau2 <- sqrt(se.Q^2/(si(1)-si(2)/si(1))^2)
  
  as <- vi * (si(1) - si(2)/si(1))
  bs <- as - (k - 1)
  se.Rb <- sqrt(se.Q^2 * (sum(as/(Q + bs)^2)/k)^2)
  se.I2 <- sqrt(se.Q^2 * (k - 1)^2/Q^4)
  se.Ri <- sqrt(se.Q^2 * (k - 1 - cv2)^2/(Q - cv2)^4)
  se.CVb <- sqrt(as.double(tau2 * (se_yr^2)/(se_yr^4) + se.tau2^2/
                             (4*y_random^2*tau2)))
  
  ## Final results
  return(list(yi = yi, vi = vi, Q = Q, tau2 = tau2, 
              Rb = Rb, I2 = I2, Ri = Ri, 
              CVb = CVb, cv2 = cv2, s2I2 = s2I2, s2Ri = s2Ri, 
              se.Q = se.Q, se.tau2 = se.tau2, se.Rb = se.Rb, 
              se.I2 = se.I2, se.Ri = se.Ri, se.CVb = se.CVb))
}

## Functions for simulating data
fun <- function(k, s2s, tau2, c){
  (sum(1/(k*s2s + tau2)) - 1/c)^2
}
simMeta <- function(R, S, RR, tau2, cv_s2s){   
  c <- tau2/(R*S)
  mean_s2s <- c*S - tau2
  var_s2s <- (mean_s2s*cv_s2s)^2
  logmu <- log(mean_s2s^2/sqrt(var_s2s + mean_s2s^2))
  logsigma  <- sqrt(log(1 + var_s2s/mean_s2s^2))
  s2s <- rlnorm(S, logmu, logsigma)  
  c2 <- optimize(fun, interval = c(0, 50000), 
                 s2s = s2s, tau2 = tau2, c = c)$minimum
  sigmas2s <- c2*s2s
  betas <- rnorm(S, log(RR), sqrt(tau2))
  beta <- rnorm(S, betas, sqrt(sigmas2s))
  return(data.frame(yi = beta, vi = sigmas2s))
}


## -----------------------------------------------------------------------------
## Simulation design


## Parameters:
RR <- 2  # Change to .5, 1, 1.5, 4
CV_beta <- c(.5, 1, 3)
tau2 <- (CV_beta*log(RR))^2
S <- c(5, 10, 20, 50, 100, 200)
cv_s2s <- c(.5, 1, 2)
R <- c(.1, .3, .5, .7, .9)
comb <- expand.grid(RR = RR, S = S, cv_s2s  = cv_s2s, R = R, tau2 = tau2)
N <- 10000

# Seed for reproducibility
set.seed(1234)

start.time <- Sys.time()
#Time difference of 33.09193 mins
for (i in seq_along(tau2)){
  
  combi <- comb[comb$tau2 == tau2[i], ]

  ## Simulating data
  listData <- lapply(split(combi, seq_along(row.names(combi))), function(x){
    lapply(as.list(1:N), function(y){
      simMeta(R = x$R, S = x$S, RR = x$RR, cv_s2s = x$cv_s2s, tau2 = x$tau2)
    })
  })
  # listData: list of length 90 of lists of length N
  
  # checkSim <- mapply(function(x, y){
  #    sapply(x, function(z) round(1/sum(1/(z$vi + y$tau2)) - with(y, tau2/(R*S)), 5))
  # }, listData, split(combi, seq_along(row.names(combi))))
  # sum(checkSim > 0)

  ## Running meta-analysis
  listModel <- lapply(listData, function(z){
    lapply(z, function(x){
      het_measures(yi = x$yi, vi = x$vi)
      # Too much time and memory
      #rma.uni(yi = yi, vi = vi, data = x, method = "DL")
    })
  })
  # listModel: list of length 90 of lists of length N (each a list of 17 elements)
  
  ## Extracting heterogeneity measures (90 lists, each with matrix Nx3)
  listMeasures <- lapply(listModel, function(z){
    do.call("rbind", lapply(z, function(x){
      unlist(x[c("Rb", "I2", "Ri", "se.Rb", "se.I2", "se.Ri")])
    }))
  })
  meanMeasure <- t(sapply(listMeasures, function(x) 
    apply(x[, c("Rb", "I2", "Ri")], 2, mean)))
  # data.frame: 90x8
  
  ## 1. Bias
  bias <- cbind(combi, apply(meanMeasure, 2, function(m) 100*(m - combi$R)/combi$R))
  biasSum <- reshape(bias, v.names = c("Rb", "I2", "Ri"), timevar = "cv_s2s",
                     idvar = c("R", "S", "tau2"), direction = "wide")

  ## 2. Coverage
  listCI <- lapply(listMeasures, function(x)
    cbind(
      x[, c("Rb", "I2", "Ri")] - qnorm(.975)*x[, c("se.Rb", "se.I2", "se.Ri")],
      x[, c("Rb", "I2", "Ri")] + qnorm(.975)*x[, c("se.Rb", "se.I2", "se.Ri")]
    ))
  coverage <- cbind(combi, do.call("rbind", mapply(function(x, r){
    100*apply(x[, 1:3] < r$R &  x[, 4:6] > r$R, 2, sum)/nrow(x)
  }, listCI, split(combi, seq_along(row.names(combi))), SIMPLIFY = FALSE)))
  coverageSum <- reshape(coverage, v.names = c("Rb", "I2", "Ri"), timevar = "cv_s2s",
                         idvar = c("R", "S", "tau2"), direction = "wide")
  
  
  ## Saving and removing objects created
  save(list = c(#"listModel", 
    "combi", "listMeasures", "meanMeasure", "biasSum", "coverageSum"),
       file = paste0("data/170728_main_sim", i, ".RData"))
  rm(list = c("bias", "biasSum", "combi", "meanMeasure", "coverage", "coverageSum", 
              "listCI", "listData", "listMeasures", "listModel"))

}
end.time <- Sys.time()
time.taken <- end.time - start.time