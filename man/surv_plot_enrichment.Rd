\name{surv_plot_enrichment}
\alias{surv_plot_enrichment}
\title{Plotting Clinical Trial Metrics for Prognostic Enrichment (Survival Outcomes)}
\description{
This function plots summaries of prognostic enrichment of clinical trials with survival outcomes, based on clinical trial metrics estimated by \code{surv_enrichment}.
}
\usage{
surv_plot_enrichment(x, km.quantiles = c(0,0.25,0.5,0.75),
                     km.range = NULL, alt.color = NULL)
}
\arguments{
  \item{x}{
    Object returned by \code{surv_enrichment}.
  }
  \item{km.quantiles}{
    Enrichment levels on which Kaplan-Meier survival estimates (Plot 1) are plotted. Defaults to four quartiles.
  }
  \item{km.range}{
    (Optional) a scalar specifying the range of time for which Kaplan-Meier survival estimates (Plot 1) are plotted. Defaults to the last time point of observation.
  }
  \item{alt.color}{
    (Optional) allows the user to specify the color of curves for clinical trial metrics (Plots 2-6). The length should match the number of trial lengths considered. Defaults to \code{ggplot2} color palette.
  }
}
\value{
A grid containing either the first 4 or 6 plots described below.
  \item{km.plot}{
    The Kaplan-Meier survival curves for specified enrichment levels. The vertical reference line(s) correspond to \code{end.of.trial} or \code{a,f}. This will be presented even if \code{method = "NNE"} was specified.
  }
  \item{prob.plot}{
    The estimated event probability (and 95\% confidence intervals) at each enrichment level.
  }
  \item{ss.plot}{
    The estimated sample size (and confidence intervals) at each enrichment level.
  }
  \item{screen.plot}{
    The estimated number of patients that need to be screened (and confidence intervals) to enroll the trial.
  }
  \item{cost.plot}{
    The estimated total cost of the trial (and confidence intervals).
  }
  \item{reduction.cost.plot}{
    The percentage of reduction in total cost comparing an enriched versus unenriched trial.
  }
  \item{summary}{
    A grid of the first 4 or all 6 plots combined together.
  }
}
\examples{
## Following the example of 'surv_enrichment':
data(SurvMarkers)
SurvMarkers$surv <- Surv(SurvMarkers$time, SurvMarkers$event)
rslt1 <- surv_enrichment(formula = surv~x2, data = SurvMarkers, hr = 0.8, end.of.trial=c(36,48),
                         cost.screening = 300, cost.keeping = c(4000,5000), cost.unit.keeping = NULL,
                         method = "KM", power = 0.9, alpha = 0.05, one.sided = F,
                         selected.biomarker.quantiles = seq(from = 0, to = 0.9, by = 0.1),
                         do.bootstrap = FALSE, print.summary.tables = FALSE)
## Truncate the range of x axis in the Kaplan-Meier plot to 0-60 months, and use colors 'salmon' and 'royalblue' for 36- and 48-month trials respectively.
plots1 <- surv_plot_enrichment(rslt1, km.quantiles = c(0,0.25,0.5,0.75),
                     km.range = 60, alt.color = c("salmon","royalblue"))
}
