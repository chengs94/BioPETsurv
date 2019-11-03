# simulating data with biomarker and survival observations
library(ggplot2)
library(survival)

sim_data <- function(n = 500, biomarker = "normal", effect.size = 0.5,
                     baseline.hazard = "constant", end.time = 10, end.survival = 0.5, shape = NULL,
                     seed = 2333){
  # effect size is log(HR) when sd(biomarker)=1
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
  X <- biom
  #colnames(X) <- "biomarker"
  #if (ncol(X)>1){
  #  colnames(X) <- c("biomarker",paste("x",seq(1:(ncol(X)-1)),sep=""))
    #b.scale <- apply(X,2,sd)
    #X.scaled <- cbind(X[,1],apply(X,2,function(x) (x-mean(x))/sd(x))[,-1])
    #beta.scaled <- beta/b.scale[-1]
  #}
  #b <- matrix(c(effect.size, beta),ncol=1)
  b <- log(effect.size)
  #Xb <- X%*%b
  Xb <- X*b
  hr <- exp(Xb)
  if (baseline.hazard=="constant"){
    lambda.surv <- -log(end.survival)/end.time
    #lambda.cens <- prob.censor*lambda.surv/(1-prob.censor)
    lambda <- lambda.surv*hr
    t.surv <- rexp(n, rate = lambda)
    #t.cens <- rexp(n, rate = lambda.cens)
    t.cens <- rep(end.time, n)
  } else{
    if (baseline.hazard=="increasing"){
      if (is.null(shape)){
        print("No Weibull shape parameter specified; defaults to shape = 2")
        k <- 2
      } else if (shape <= 1){
        stop("Weibull shape should >1 for an increasing baseline hazard")
      }
      else k <- shape
    }
    if (baseline.hazard=="decreasing"){
      if (is.null(shape)){
        k <- 0.5
        print("No Weibull shape parameter specified; defaults to shape = 0.5")
      } else if (shape >= 1){
        stop("Weibull shape should <1 for a decreasing baseline hazard")
      }
      else k <- shape
    }
    b.bl <- end.time*(-log(end.survival))^(-1/k)
    lambda.surv <- b.bl^(-k)
    lambda <- lambda.surv*hr
    b <- lambda^(-1/k)
    t.surv <- rweibull(n, shape = k, scale = b)
    t.cens <- rep(end.time, n)
  }
  t.obs <- apply(cbind(t.surv,t.cens), 1, min)
  event <- apply(cbind(t.surv,t.cens), 1, function(vec) as.numeric(vec[1]<vec[2]))
  t.obs[t.obs>end.time] <- end.time
  event[t.obs==end.time] <- 0
  #dat <- cbind(X,t.obs,event,t.surv,t.cens)
  dat <- cbind(X,t.obs,event,t.surv)
  dat <- as.data.frame(dat)
  colnames(dat) <- c("biomarker","time.observed","event","time.event")
  #colnames(dat)[(ncol(X)+1):ncol(dat)] <- c("time.observed","event","time.event","time.censoring")

  # plotting the K-M curve
  cols <- gray.colors(7)
  km.quantiles <- c(0, 0.25, 0.5, 0.75)
  km.all <- survfit(Surv(dat$time.observed, dat$event)~1, error="greenwood")
  dat1 <- as.data.frame(seq(0,max(km.all$time),by=max(km.all$time)/500))
  colnames(dat1) <- "time"
  for (j in 1:length(km.quantiles)){
    q <- quantile(dat$biomarker,prob=km.quantiles[j])
    sobj <- Surv(dat$time.observed, dat$event)[dat$biomarker>=q]
    km <- survfit(sobj~1,error="greenwood")
    survfun <- stepfun(km$time, c(1, km$surv))
    dat1 <- cbind(dat1, survfun(dat1[,1]))
    colnames(dat1)[j+1] <- paste(j,"surv",sep=".")
  }
  dat1 <- reshape(dat1, direction = 'long', timevar = 'level.enrichment',
                 varying=list(grep("surv", colnames(dat1), value=T)),
                 times = as.character(km.quantiles),
                 v.names = c("surv"),
                 idvar='time')
  g <- ggplot(dat1,aes(x=time, y=surv, colour=level.enrichment)) +
    geom_line(size=1) + ylim(0,1) +
    labs(title ="Kaplan-Meier survival curves",
         x = "time", y = "survival estimate", color = "enrichment level") +
    scale_color_manual(labels = as.character(km.quantiles), values = cols[1:4]) +
    theme(plot.title = element_text(hjust = 0.5), legend.position="bottom")

  plot(g)
  return(list(data = dat, km.plot = g))
}

