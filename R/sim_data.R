# simulating data with biomarker and survival observations

sim_data <- function(n = 500, covariates = NULL, beta = NULL, biomarker = "normal", effect.size = 0.25,
                     baseline.hazard = "constant", end.time = 10, end.survival = 0.5, prob.censor = 0.2,
                     seed = 2333){
  # effect size is log(HR) when sd(biomarker)=1
  if (!is.null(covariates) & !is.null(beta)){
    if (nrow(covariates)!=n | ncol(covariates)!=length(beta))
      stop("Invalid dimension of covariates (# rows should match the sample size n, # cols should match the dimension of beta)")
  }
  if (!baseline.hazard %in% c("constant","increasing","decreasing")){
    stop("Invalid type of baseline hazard (should be constant/increasing/decreasing)")
  }
  if (!biomarker %in% c("normal","lognormal")){
    stop("Invalid distribution of biomarker (should be normal/lognormal)")
  }
  if (end.survival <=0 & end.survival >=1){
    stop("end.survival should be between 0 and 1")
  }

  ##### simulating the data ######
  set.seed(seed)
  biom <- rnorm(n)
  if (biomarker=="lognormal") biom <- (exp(biom)-mean(exp(biom)))/sd(exp(biom))
  X <- cbind(biom, covariates)
  colnames(X)[1] <- "biomarker"
  if (ncol(X)>1){
    colnames(X) <- c("biomarker",paste("x",seq(1:(ncol(X)-1)),sep=""))
    #b.scale <- apply(X,2,sd)
    #X.scaled <- cbind(X[,1],apply(X,2,function(x) (x-mean(x))/sd(x))[,-1])
    #beta.scaled <- beta/b.scale[-1]
  }
  b <- matrix(c(effect.size, beta),ncol=1)
  Xb <- X%*%b
  if (baseline.hazard=="constant"){
    lambda.surv <- -log(end.survival)/end.time
    lambda.cens <- prob.censor*lambda.surv/(1-prob.censor)
    #lambda.cens <- (1-prob.censor)*lambda.surv/prob.censor
    hr <- exp(Xb)
    lambda <- lambda.surv*hr
    t.surv <- rexp(n, rate = lambda)
    t.cens <- rexp(n, rate = lambda.cens)
    #t.cens[t.cens>end.time] <- end.time
    t.obs <- apply(cbind(t.surv,t.cens), 1, min)
    event <- apply(cbind(t.surv,t.cens), 1, function(vec) as.numeric(vec[1]<vec[2]))
    t.obs[t.obs>end.time] <- end.time
    event[t.obs==end.time] <- 0
    dat <- cbind(X,t.obs,event,t.surv,t.cens)
    dat <- as.data.frame(dat)
    colnames(dat)[(ncol(X)+1):ncol(dat)] <- c("time.observed","event","time.event","time.censoring")
  }
  return(dat)
}


#temp=sim_data()
#mean(temp$event)
#mean(temp$time.event<temp$time.censoring)
