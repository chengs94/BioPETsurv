# BioPETsurv: Biomarker Prognostic Enrichment Tool for clinical trials with survival outcomes

Prognostic Enrichment is a clinical trial strategy of evaluating an intervention in a patient population with a higher rate of the unwanted event than the broader patient population (R. Temple (2010) <DOI:10.1038/clpt.2010.233>). A higher event rate translates to a lower sample size for the clinical trial, which can have both practical and ethical advantages. The package ```BioPET-Surv``` provides tools to evaluate biomarkers for prognostic enrichment of clinical trials with survival or time-to-event outcomes. An associated R shiny webtool (with simplified functionality) can be found [here](https://chengs94.shinyapps.io/biopetsurv/).

Key functions of this package are:

* ```sim_data```: Simulate a dataset containing biomarker and survival observations
* ```surv_enrichment```: Estimate trial characteristics at different levels of enrichment, given a biomarker (which can be a single biomarker or a composite)
* ```surv_plot_enrichment```: Visualize trial characteristics returned by ```surv_enrichment```


**Update history**

v4, 9/19/2019

Added the functionality of simulating datasets containing biomarker and survival observations. The R function currently allows for constant baseline hazard.

v3, 5/21/2019:

Incorporated an alternative method for calculating event rates. The method comes from [Heagerty et al (2000)](https://www.ncbi.nlm.nih.gov/pubmed/10877287), and uses a kernel smoothed version of Kaplan-Meier survival estimators. This method allows the censoring process to be dependent on the biomarker, and guarantees monotone ROC curves (if the user is interested in the prognostic capacity of a biomarker represented by time-dependent ROC curves).
