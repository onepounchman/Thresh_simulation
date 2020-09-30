library(simcausal)
library(dplyr)

source("simu_functions.R")

# create DAG object
setD=generate_covariates()

# Y follows a logit model
dat_logit=simu_logit_cov(0.1,n=1000,p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09),p0=0.2,p1=1,seed=1,setD=setD)
mean(dat_logit$simmed_data$Y)

# Y follows a probit model
dat_probit=simu_probit_cov(0.1,n=1000,p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09),p0=0.2,p1=1,seed=1,setD=setD)
mean(dat_probit$simmed_data$Y)

# Y follows a step model
dat_step=simu_step(0.1,n=1000,p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09),p0=0.2,p1=1,seed=1,setD)
mean(dat_step$simmed_data$Y)

# Y follows a scaled logit model
dat_scale=simu_logit_cov(0.1,lambda=0.5,n=1000,p=c(0.009, 0.01, 0.03, 0.05, 0.07, 0.09),p0=0.2,p1=1,seed=1,setD=setD)
mean(dat_scale$simmed_data$Y)

