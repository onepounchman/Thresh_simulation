library(resample)
library(simcausal)
library(dplyr)
library(tmleThresh)
library(tmle3)
library(sl3)
library(tibble)
library(ggplot2)

source('simu_functions_v3.R')


#logit_dat_unbiased<-simu(Risk=0.3,p0=1,p1=1,cov='non_confound',model='logit') 

# risk level 0.3: mean(data$Y) close to 0.3

p0=1
p1=1

D <- DAG.empty()
D <- D +
  node("W1", distr = "runif", min = -1, max = 1) +
  node("W2", distr = "rnorm", mean = 0) +
  node("W3full", distr = "rexp", rate = 1) +
  node("W3", distr = "rconst", const = min(W3full, 3)) +
  node("Sfull", distr = "rgamma", shape = 2.5) +
  node("S", distr = "rconst", const=min(Sfull, 13)/13) +
  node("PY", distr = "rconst", const = logit(0.2*(-2.2-7*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1)))) +
  node("Y", distr = "rbinom",size =1, prob = PY)+
  node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
D <- set.DAG(D)

# true value
true<-sim_thre(D,30000,rndseed=1,p=0.5,cutoffs=seq(0.1,0.9,by=0.1))$true

unbiased_2000<-sim_sample(logit_dat_unbiased,n = 2000,rndseed=0)

multiple<-tmle_esti(unbiased_2000,seq(0.1,1,by=0.1))
single<-tmle_esti(unbiased_2000,0.5)
dov<-Donov_esti(unbiased_2000,0.5,unbiased_2000$weights)


## create tables
# single cut offs
tabl1<-get_table_single(D,n=2000,times=5,cutoffs=seq(0.1,0.5,by=0.1))

# multiple cut offs
tabl2<-get_table(D,n=2000,times=5,cutoffs=seq(0.1,0.5,by=0.1))




