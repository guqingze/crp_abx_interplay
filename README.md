# Interplay between CRP responses and antibiotic prescribing

R code for paper "Interplay between CRP responses and antibiotic prescribing in suspected bloodstream infections".

## Input data
The data analysed are available from the Infections in Oxfordshire Research Database (https://oxfordbrc.nihr.ac.uk/research-themes/modernising-medical-microbiology-and-big-infection-diagnostics/iord-about/), subject to an application and research proposal meeting on the ethical and governance requirements of the Database.

## Software and packages
Analyses were performed using statistical software R, version 4.3.2 (R Project for Statistical Computing). 
Data management used the ‘data.table’ package (version 1.14.0). 
Model fitting was performed using the ‘nnet’ package (version 7.3-19), the ‘lme4’ package (version 1.1-27.1, using restricted maximum likelihood estimation (REML) with the ‘t-tests use Satterthwaite’s method’ approach for approximating degrees of freedom, as implemented in the ‘lmerModLmerTest’ function), the ‘stats’ package (version 4.3.2), and the ‘splines’ package (version 4.3.2).
Graphs were produced using the ‘ggplot2’ package (version 3.5.1), the ‘ggsurvfit’ package (version 1.1.0), the ‘patchwork’ package (version 1.2.0), and the ‘viridis’ package (version 0.6.1).

## Code
- computeRxDecision.R: R code for computing antibiotic prescribing decisions based on differences in this pre-defined antibiotic ranking, number of different antibiotics administered, and route of administration.
- multinom_abx_after_crp.R: R code for multinomial logistic regression analysis of antibiotic prescribing after CRP centile changes.
- lmer_crp_after_abx.R: R code for linear mixed-effects regression analysis of CRP centile changes after antibiotic prescribing changes.
- mortality_early_crp.R: R code for analysis of the impact of early CRP centile changes on 5-30 day all-cause mortality.
- mortality_auc_compare.R: R code for comparison of the area under the receiver operating characteristic curve (AUC) for various models for 5-30 day all-cause mortality, including covariates only, centile changes, absolute CRP change, and log CRP change per day (days 1-4), and combinations of covariates with these.