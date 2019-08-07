\name{surv_enrichment}
\alias{surv_enrichment}
\title{Prognostic Enrichment for Survival Outcomes}
\description{
This function evaluates biomarkers for prognostic enrichment of clinical trials with survival outcomes, using data with biomarker values and survival observations.
}
\usage{
surv_enrichment(formula, data, hr = 0.8, end.of.trial=NULL, a=NULL, f=NULL,
     method = "KM", lambda = 0.05,
     cost.screening = NULL, cost.keeping = NULL, cost.unit.keeping = NULL,
     power = 0.9, alpha = 0.05, one.sided = F,
     selected.biomarker.quantiles = seq(from = 0, to = 0.95, by = 0.05),
     do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333,
     print.summary.tables = FALSE)
}
\arguments{
  \item{formula}{
    Object of class "formula", in the form "outcome ~ predictors", where the outcome is a survival object as returned by function "Surv()". The predictor(s) include the biomarker(s) of interest and/or other covariates. If multiple predictors are inputed, this function will consider a composite biomarker based on a Cox model with all predictors.
  }
  \item{data}{
    Data frame containing the survival outcome (as returned by "Surv()") and predictors specified in the "formula" argument. Observations with missing value(s) will be dropped.
  }
  \item{hr}{
    The hazard ratio (comparing the treatment versus control group) that the trial seeks to detect. Should be  a number between 0 and 1.
  }
  \item{end.of.trial}{
    A scalar or a vector that gives the duration(s) of fixed-length trial(s) that are being considered. Either "end.of.trial" or a combination of "a" and "f" must be specified.
  }
  \item{a}{
    A scalar specifying the accrual period of the trial. It is assumed that patients are recruited at a uniform rate during this period.
  }
  \item{f}{
    A scalar specifying the follow-up period of the trial, that is, all patients are followed during this period right after the accrual period.
  }
  \item{method}{
    Character specifying which
  }
}
\examples{
hello()
}
