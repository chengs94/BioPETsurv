# demo with assess: new

load("C:/Users/Si Cheng/OneDrive - UW/19winter/biomarker/assess_aki/SurvPET_data.RData")
load("C:/Users/Si Cheng/OneDrive - UW/19winter/biomarker/SurvPET/sim_data.RData")
#source("C:/Users/Si Cheng/OneDrive - UW/19winter/biomarker/SurvPET/R/funcs.R")
source("C:/Users/Si Cheng/OneDrive - UW/19winter/biomarker/SurvPET/R/surv_enrichment.R")

## C-stat using simulated dataset
library(survcomp)
cstat1=concordance.index(sim.data$x1, sim.data$surv[,1], sim.data$surv[,2])$c.index
cstat2=concordance.index(sim.data$x2, sim.data$surv[,1], sim.data$surv[,2])$c.index

rslt.short.trial <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=2.8, a=NULL, f=NULL,
                                    cost.screening = 50, cost.keeping = 100, cost.unit.keeping = NULL,
                                    power = 0.9, alpha = 0.05, one.sided = F,
                                    selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                                    do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333)


# using the original dataset
# no accrual+follow up; uniform cost for trial; two lengths of trial
rslt1 <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=c(36,48), a=NULL, f=NULL,
                         cost.screening = 50, cost.keeping = c(800,1000), cost.unit.keeping = NULL,
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                         do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333)
plots1 <- surv_plot_enrichment(rslt1)

# same as above; use RMST to compute cost
rslt1.5 <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=c(36,48), a=NULL, f=NULL,
                         cost.screening = 50, cost.keeping = NULL, cost.unit.keeping = c(25,25),
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                         do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333)
plots1.5 <- surv_plot_enrichment(rslt1.5)


# no accrual+follow up; cost by unit time; alternative color
rslt2 <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=c(36,48), a=NULL, f=NULL,
                         cost.screening = 50, cost.keeping = NULL, cost.unit.keeping = 25,
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                         do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333)
plots2 <- surv_plot_enrichment(rslt2, alt.color = c("salmon","royalblue"))

# accrual + follow up; uniform cost, no bootstrap sd
rslt3 <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=NULL, a=24, f=24,
                         #cost.screening = 50, cost.keeping = 900, cost.unit.keeping = NULL,
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                         do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333)
plots3 <- surv_plot_enrichment(rslt3)

# accrual + follow up; varying cost, with bootstrap sd
start <- Sys.time()
rslt4 <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=NULL, a=24, f=24,
                         cost.screening = 50, cost.keeping = NULL, cost.unit.keeping = 25,
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                         do.bootstrap = T, n.bootstrap = 500, seed = 2333)
Sys.time()-start # 4.86min
plots4 <- surv_plot_enrichment(rslt4)

# combined biomarker; accrual + follow up; uniform cost; bootstrap sd
start <- Sys.time()
rslt5 <- surv_enrichment(formula = survobj~V3M_YKL40+V3M_IL18, data, hr = 0.8, end.of.trial=NULL, a=24, f=24,
                         cost.screening = 50, cost.keeping = 900, cost.unit.keeping = NULL,
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.98, by = 0.02),
                         do.bootstrap = T, n.bootstrap = 500, seed = 2333)
Sys.time()-start
# specified range for K-M plot
plots5 <- surv_plot_enrichment(rslt5, km.range = 60)

# using NNE estimate; compare with rslt1
rslt6 <- surv_enrichment(formula = survobj~V3M_YKL40, data, hr = 0.8, end.of.trial=48, a=NULL, f=NULL,
                         method="NNE", lambda=0.05,
                         cost.screening = 50, cost.keeping = 1000, cost.unit.keeping = NULL,
                         power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.9, by = 0.1),
                         do.bootstrap = FALSE)
plots6 <- surv_plot_enrichment(rslt6)


# simulated data based on ASSESS
load("C:/Users/Si Cheng/OneDrive - UW/19winter/biomarker/assess_aki/SurvPET_data.RData")
set.seed(2333)
sim.data <- na.omit(data)
sim.data$x1 <- sim.data$V3M_YKL40 + rnorm(nrow(sim.data),mean=0,sd=15)
sim.data$t <- sim.data$time_ckdipe + rnorm(nrow(sim.data),mean=0,sd=3)
sim.data$t[sim.data$t<0] <- 1
sim.data$newsurv <- Surv(sim.data$t, as.numeric(sim.data$ckdipe))
sim.data$x2 <- -sim.data$time_ckdipe + rnorm(nrow(sim.data),mean=0,sd=20)
sim.data <- sim.data[,c(13,16,15)]
sim.data$x1 <- (sim.data$x1-mean(sim.data$x1))/sd(sim.data$x1)
sim.data$x2 <- (sim.data$x2-mean(sim.data$x2))/sd(sim.data$x2)
colnames(sim.data) <- c("x1","x2","surv")
