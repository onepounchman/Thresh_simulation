

simu<-function(Risk,p0,p1,cov,model){
# Arguments:
  # Risk: disease risk
  # n: sample size
  # p0: immune response sampling rate for control group case-control study design
  # p1: immune response sampling rate for case group case-control study design
  # cov: whether inclduing covariates and how to use covariates
  # model: for P(Y|S=v,W=w)
#------------------------------------------------------------------------
# output: A DAG subject
#----------------------------------
  D<-generate_covariates(cov)
  D<-generate_S(D,cov)
  if(cov=='non_covariate'){
    fun_y<-eval(parse(text='generate_Y_noncov'))
    D<-fun_y(Risk,D,model,p0,p1)
  }
  if(cov=='non_confound'){
    fun_y<-eval(parse(text='generate_Y_nonconfound'))
    D<-fun_y(Risk,D,model,p0,p1)
  }
  if(cov=='confound'){
    fun_y<-eval(parse(text='generate_Y_confound'))
    D<-fun_y(Risk,D,model,p0,p1)
  }
  D
}



rmults <- function(n, size, prob) {
  size = size[1]
  prob = prob[1]
  apply(rmultinom(n = n, 1:size, rep(prob, size))==1, 2, which)

}

logit <- function(x){1/(1+exp(-x))}
scale <- function(x){0.7*1/(1+exp(-x))}
probit<-function(x){pnorm(x)}


  
generate_covariates<-function(cov){
    # create a DAG object
  D <- DAG.empty()
  if(cov=='non_covariate'){
    D <- D +
      node("A", distr = "rmults", size = 10 , prob = 1/10) +
      node("B", distr = "rmults", size = 10 , prob = 1/10) +
      node("C", distr = "rmults", size = 10 , prob = 1/10) +
      node("W1", distr = "rconst", const = A -3) +
      node("W2", distr = "rconst", const = B -3) +
      node("W3", distr = "rconst", const = C -3) 
  }
  else{
    D <- D +
      node("W1", distr = "runif", min = -1, max = 1) +
      node("W2", distr = "rnorm", mean = 0) +
      node("W3full", distr = "rexp", rate = 1) +
      node("W3", distr = "rconst", const = min(W3full,3))
  }
  D
}




generate_S<-function(D,cov){
  if(cov=='non_covariate') {D<-D+node("S", distr = "rgamma",  shape=4,rate=1)}
  else if(cov=='non_confound') {
    D<-D + node("Sfull", distr = "rgamma", shape = abs(2.5)  ) +
      node("S", distr = "rconst", const = min(Sfull, 13)/13)}
  else if(cov=='confound') {
    D<-D+
      node("Sfull", distr = "rgamma", shape = abs(2 + W3full + W3full*W2 + W2*exp(W1) + cos(W2) + sin(W1)))  +
      node("S", distr = "rconst", const = min(Sfull, 13)/13)
    } 
  D
}

generate_Y_nonconfound<-function(Risk,D,model,p0,p1){
      fun=eval(parse(text=model))
      beta0<-seq(-10,10,0.1)
      if(model=='step') beta0<-seq(0.01,0.99,0.01)
      risk<-rep(0,length(beta0))
      for(i in 1:length(beta0)){
        new_D<-D
        if(model!='step'){
          #1/4*(eval(beta0[i])-S/2-W1*(W2 >0) - W2*(W1>0) + W3*S - sin(W2) + cos(W1) logit
          #
          new_D<-D+node("PY", distr = "rconst", const = fun(0.2*(eval(beta0[i])-7*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1)))) +
            node("Y", distr = "rbinom", size =1, prob = PY)
          }
        else{
          new_D<-D+node("PY",distr = "rconst", const = eval(beta0[i]))+
            node("Y", distr = "rbinom", size =1, prob = PY)
        }
        setD<-set.DAG(new_D)
        data<- sim(setD, n = 50000,rndseed=1)
        risk[i]<-mean(data$Y)
      }
      mydiff <- abs(risk-Risk)
      mymin <- min(mydiff)
      select <- head(beta0[mydiff==mymin],1)
      if(model!='step'){
        D<-D+node("PY", distr = "rconst", const = fun(0.2*(eval(select)-7*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1))))+
          node("Y", distr = "rbinom", size =1, prob = PY)
        }
      else{
        D<-D+node("PY",distr = "rconst", const = eval(select))+
          node("Y", distr = "rbinom", size =1, prob = PY)
        }
      D<-D+node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
      set.DAG(D)
}

generate_Y_noncov<-function(Risk,D,model,p0,p1){
  fun=eval(parse(text=model))
  beta0<-seq(-10,10,0.1)
  if(model=='step') beta0<-seq(0.01,0.99,0.01)
  risk<-rep(0,length(beta0))
  for(i in 1:length(beta0)){
    new_D<-D
    if(model!='step'){
      new_D<-D+node("PY", distr = "rconst", const = fun(eval(beta0[i])-5*S))+
        node("Y",distr='rbern',prob=PY)
      }
    else{
      new_D<-D+node("PY",distr = "rconst", const = eval(beta0[i]))+
        node("Y",distr='rbern',prob=PY)
      }
    setD<-set.DAG(new_D)
    data<- sim(setD, n = 50000,rndseed=1)
    risk[i]<-mean(data$Y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  if(model!='step'){
    D<-D+node("PY", distr = "rconst", const = fun(eval(select)-5*S))+
      node("Y",distr='rbern',prob=PY)
    }
  else{D<-D+node("PY",distr = "rconst", const = eval(select))+
    node("Y",distr='rbern',prob=PY)
  }
  
  D<-D+node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
  setD<-set.DAG(D)
}

generate_Y_confound<-function(Risk,D,model,p0,p1){
  fun=eval(parse(text=model))
  beta0<-seq(-10,10,0.1)
  if(model=='step') beta0<-seq(0.01,0.99,0.01)
  risk<-rep(0,length(beta0))
  for(i in 1:length(beta0)){
    new_D<-D
    if(model!='step'){
      new_D<-D+node("PY", distr = "rconst", const = fun(0.2*(eval(beta0[i])-4*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1))))+
        node("Y", distr = "rbinom", size =1, prob = PY)
    }
    else{
      new_D<-D+node("PY",distr = "rconst", const = eval(beta0[i]))+
        node("Y",distr='rbern',prob=PY)
    }
    setD<-set.DAG(new_D)
    data<- sim(setD, n = 50000,rndseed=1)
    risk[i]<-mean(data$Y)
  }
  mydiff <- abs(risk-Risk)
  mymin <- min(mydiff)
  select <- head(beta0[mydiff==mymin],1)
  if(model!='step'){
    D<-D+node("PY", distr = "rconst", const = fun(0.2*(eval(select)-4*S + W1+W2+2*W1*(W2 >0) - 2*W2*(W1>0) + W3 * (S-0.1) - W3*cos(W1))))+
      node("Y", distr = "rbinom", size =1, prob = PY)
  }
  else{
    D<-D+node("PY",distr = "rconst", const = eval(select))+
      node("Y",distr='rbern',prob=PY)
  }
  D<-D+node('G',distr='rbern',prob=(Y==1)*eval(p1)+(Y==0)*eval(p0))
  set.DAG(D)
}


sim_sample<-function(object,n,rndseed=1){
data<-sim(object,n=n ,rndseed=rndseed)

data<- data %>% mutate(S=ifelse(G==1, S, NA))

p0_hat <- mean(!is.na(data[data$Y==0,]$S))
p1_hat <- mean(!is.na(data[data$Y==1,]$S))

IPW <- function(x){
  term <- ifelse(x==0, p0_hat, p1_hat)
  1/term
}

weights <- unlist(lapply(data$Y, IPW))

data$weights=weights

output<-data

return(output)
}




sim_thre<-function(object,n,rndseed=1,p,cutoffs){
  data<-sim(object,n=n ,rndseed=rndseed)
  
  truethres <- rep(NA, length=length(p))
  v <- seq(0,8,by=.001)
  nv <- length(v)
  risks <- vector(, length=nv)
  for(i in 1:nv){
    risks[i]<-mean(data[data$S>=v[i],]$PY)
  }
  for(j in 1:length(p)){
    mydiff <- abs(risks-p[j])
    mymin <- min(mydiff,na.rm = TRUE)
    index <- which(mydiff==mymin)
    truethres[j] <- min(v[index])
  }
  print(data)
  true<-sapply(cutoffs,function(x) mean(data[data$S>=x,]$PY))

  data$G=1
  data<- data %>% mutate(S=ifelse(G==1, S, NA))
  #p0_hat <- mean(!is.na(data[data$Y==0,]$S))
  #p1_hat <- mean(!is.na(data[data$Y==1,]$S))

  
  output<-list(data, truethres,true=true)

  names(output)[1] <- "simmed_data"
  names(output)[2] <- "true_thresholds"
  names(output[[2]]) <- as.character(p)
  output
}

IPW <- function(x){
  term <- ifelse(x==0, p0_hat, p1_hat)
  1/term
}


## estimation functions

#---------------------------------------------------------------
## functions
get_table<-function(object,n,times,cutoffs=seq(0.1,1,by=0.1),p=c(0.01,0.03,0.05,0.07,0.09)){
  true=sim_thre(object,30000,rndseed=1,p,cutoffs)$true
  table<-0
  for(i in 1:times){
    sample<-sim_sample(object,n =n,rndseed=i)
    sample_complete<-sample[complete.cases(sample),]
    donov<-Donov_esti(sample_complete,cutoffs,sample_complete$weights)
    tmle<-tmle_esti(sample_complete,cutoffs)
    table<-table+cbind(tmle,donov)
  }
  table<-table/times
  table<-add_column(table, true, .after = 1)
  round(table,3)
}

get_table_single<-function(object,n,times,cutoffs=seq(0.1,1,by=0.1),p=c(0.01,0.03,0.05,0.07,0.09)){
  true=sim_thre(object,30000,rndseed=1,p,cutoffs)$true
  table<-0
  for(i in 1:times){
    sample<-sim_sample(object,n =n,rndseed=i)
    sample_complete<-sample[complete.cases(sample),]
    list<-sapply(cutoffs,function(x){
      donov<-Donov_esti(sample_complete,x,sample_complete$weights)
      tmle<-tmle_esti(sample_complete,x)
      single<-cbind(tmle,donov)
    },simplify = FALSE)
    table_single<-Reduce(rbind,list)
    table<-table+table_single
  }
  table<-table/times
  table<-add_column(table, true, .after = 1)
  round(table,3)
}



Donov_esti<-function(data,thresholds,weights){
  get_estimates <- function(data, marker_var, outcome_var, thresholds, weights = rep(1, nrow(data))) {
    weights <- weights/sum(weights)
    data <- as.data.table(data)
    IC_list <- list()
    est_list <- list()
    for(thresh in thresholds) {
      meets_thresh <- as.numeric(data[[marker_var]] >= thresh) 
      cdf <- weighted.mean(meets_thresh, weights)
      EY <- weighted.mean(data[[outcome_var]] * meets_thresh, weights) / cdf
      thresh <- as.character(thresh)
      est_list[[thresh]] <- EY
      IC_list[[thresh]] <- meets_thresh/cdf * (data[[outcome_var]] - EY)
      
    }
    return(list(est = unlist(est_list),
                IC = do.call(cbind, IC_list) ))
  }
  #thresholds <- quantile(data[["S"]], seq(0.1, 0.9, length.out = 10))
  
  out <- get_estimates(data, "S", "Y", thresholds,weights)
  est <- out$est
  IC <- out$IC
  radius <- 1.96 * sqrt(resample::colVars(IC))/sqrt(nrow(data))
  lower <- est - radius
  upper <- est + radius
  table<-data.frame(DOV_est=est,DOV_lower=lower, DOV_upper=upper)
  rownames(table)<-seq(length(thresholds))
  table
}



tmle_esti<-function(data,thresholds){
  
  threshold_function = function(A) {thresholds}
  library(tmle3)
  #lrnr <- Lrnr_pooled_hazards$new(Lrnr_xgboost$new())
  lrnr_Y <- Lrnr_xgboost$new()
  lrnr_A <- Lrnr_xgboost$new()
  # lrnr_Y <-Lrnr_polspline$new()
  # lrnr_A <-Lrnr_polspline$new()
  tmle_spec <- tmle3_Spec_Threshold$new(method = "Psi_W",threshold_function=threshold_function)
  learner_list <- list("A" = lrnr_A, "Y" =lrnr_Y)
  node_list <- list("W" = c('W1','W2','W3'), "A" = "S", "Y" = "Y", weights = "weights")
  start_time <- proc.time()
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)
  task_time <- proc.time()
  initial_likelihood <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)
  likelihood_time <- proc.time()
  updater <- tmle_spec$make_updater(maxit = 100, verbose = T)
  targeted_likelihood <- tmle_spec$make_targeted_likelihood(initial_likelihood, updater)
  tmle_params <- tmle_spec$make_params(tmle_task, targeted_likelihood)
  
  suppressWarnings(targeted_likelihood$updater$update_step(targeted_likelihood, tmle_task))
  estimates <- tmle_params[[1]]$estimates(tmle_task)
  # This should be very small
  #colMeans(estimates$IC)
  # Estimates of Psi
  #estimates$psi
  # Estimates and confidence bounds
  #summary_from_estimates(tmle_task, list(estimates))
  
  sumry <- summary_from_estimates(tmle_task, list(estimates))
  #sumry
  lower <- sumry$lower
  upper <- sumry$upper
  est <- sumry$tmle_est
  plot_data <- data.frame(cutoffs =thresholds,est = est, lower = lower, upper = upper)
  cutoffs = round(targeted_likelihood$factor_list$Y$learner$cutoffs,3)
  #list(data=plot_data,wei=matrix(targeted_likelihood$get_likelihood(tmle_task, "A"), ncol = 10) )
  plot_data
} 
