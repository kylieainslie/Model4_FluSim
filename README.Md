---
title: "Model4_FluSim"
subtitle: "A Dynamic Model for Evaluation of the Bias of Influenza Vaccine Effectiveness Estimates from Observational Studies"
author: "Kylie Ainslie"
date: "2021-02-01"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Model4_FluSim}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Introduction

## Model Description
We present a dynamic model consisting of five steps. Below we define the model steps, the associated variables (Table 2), and the probabilities determining each variable’s distribution (Table 3). All variables are defined for each member of the study population, and we allow some variables to change over time (measured in weeks). Figure 1 illustrates the possible sources of confounding and bias present in studies designed to evaluate influenza VE [15,16]. Model assumptions are shown in Table 1.

**Step 1:** Covariates. We assume that people within the population can be classified with a health status (X) of either “healthy” or “frail” and a health awareness (U) of either “high” or “low”.

**Step 2:** Vaccination. We consider the vaccination scenario where an individual is considered vaccinated if they received the vaccine at least 14 days prior to the study onset (V = 1) or remains unvaccinated throughout the study (V = 0).

**Step 3:** Influenza and non-influenza ARI. During the influenza season, a person may become infected with an influenza virus and develop influenza ARI. Regardless of influenza infection, a person may develop one or more non-influenza ARIs. We define a variable Yj for the illness/infection status in week j as follows: Yj = 0 for no ARI, Yj = 1 for non-influenza ARI, and Yj = 2 for influenza ARI. If a person has both non-influenza ARI and influenza ARI in the same week, we consider them as influenza ARI (i.e., Yj = 2). The distribution of Yj may depend on the person’s vaccination (V) and health (X) status.

**Step 4:** Seeking medical care for ARI. A person with an ARI in week j may seek medical care (Mj). The probability of seeking medical care depends on Yj, as only those individuals who have an ARI may seek medical care, and it may be different for influenza ARI and non-influenza ARI patients. This probability may also depend on V and U.

**Step 5:** Testing for influenza infection. We assume that each person who seeks medical care for ARI is tested for influenza infection. Let Tj denote the binary test result, where T_j=1 or 0 for positive or negative, respectively.

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Table 1. Model assumptions</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Assumption </th>
   <th style="text-align:left;"> Description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> a </td>
   <td style="text-align:left;"> Each person is assigned two binary covariates: health status and health awareness, which are constant over time </td>
  </tr>
  <tr>
   <td style="text-align:left;"> b </td>
   <td style="text-align:left;"> Everyone is either effectively vaccinated prior to the study or remains unvaccinated throughout the study, and vaccination status is determined without error </td>
  </tr>
  <tr>
   <td style="text-align:left;"> c </td>
   <td style="text-align:left;"> Everyone in the study population is susceptible at the beginning of the study </td>
  </tr>
  <tr>
   <td style="text-align:left;"> d </td>
   <td style="text-align:left;"> A person can only have one influenza ARI during the study </td>
  </tr>
  <tr>
   <td style="text-align:left;"> e </td>
   <td style="text-align:left;"> A person can have at most one non-influenza ARI per week </td>
  </tr>
  <tr>
   <td style="text-align:left;"> f </td>
   <td style="text-align:left;"> The probabilities of influenza ARI and non-influenza ARI do not depend on a person’s health awareness given his/her health status </td>
  </tr>
  <tr>
   <td style="text-align:left;"> g </td>
   <td style="text-align:left;"> The probability of seeking medical care does not depend on a person’s health status given his/her health awareness </td>
  </tr>
  <tr>
   <td style="text-align:left;"> h </td>
   <td style="text-align:left;"> Influenza test specificity and sensitivity do not depend on health status, health awareness, or vaccination status, given a person’s influenza infection status </td>
  </tr>
  <tr>
   <td style="text-align:left;"> i </td>
   <td style="text-align:left;"> Every person who seeks medical care is tested for influenza infection </td>
  </tr>
  <tr>
   <td style="text-align:left;"> j </td>
   <td style="text-align:left;"> A person with no ARI does not seek medical care </td>
  </tr>
</tbody>
</table>

## Simulation Program

### Input File
<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Table 2. Input parameters</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Title </th>
   <th style="text-align:right;"> Value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Number of simulations </td>
   <td style="text-align:right;"> 100.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seed </td>
   <td style="text-align:right;"> 1234.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> All-or-none vaccine (Leaky vaccine if 'no') </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Output files (1=yes or 0=no for each) </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Incidence-weekly each simulation </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Incidence-season each simulation </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Incidence-weekly overall </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Incidence-season overall </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Outcomes file </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Add timestamp to output file names </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Year of beginning of study </td>
   <td style="text-align:right;"> 2010.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Month of beginning of study </td>
   <td style="text-align:right;"> 12.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of weeks in the study </td>
   <td style="text-align:right;"> 5.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of strata </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sizes of strata </td>
   <td style="text-align:right;"> 1000.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of X=1, U=1 (pi_{11}) </td>
   <td style="text-align:right;"> 0.400 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of X=1, U=0 (pi_{10}) </td>
   <td style="text-align:right;"> 0.400 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of X=0, U=1 (pi_{01}) </td>
   <td style="text-align:right;"> 0.100 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of X=0, U=0 (pi_{00}) </td>
   <td style="text-align:right;"> 0.100 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vaccination coverage for X=1, U=1 (alpha_{11}) </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vaccination coverage for X=1, U=0 (alpha_{10}) </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vaccination coverage for X=0, U=1  (alpha_{01}) </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vaccination coverage for X=0, U=0  (alpha_{00}) </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 0.200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Initial number of infected persons </td>
   <td style="text-align:right;"> 30.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probabilities of FARI to unprotected for X=1 (gamma_01) by week </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.008 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 0.010 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0.015 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 0.016 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FARI probabilities multipliers for X=0 (phi_gamma) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FARI probabilities multipliers for vaccinated or protected persons (theta_gamma) </td>
   <td style="text-align:right;"> 0.563 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probabilities of onset of NFARI for X=1 in unvaccinated or unprotected persons (beta_01) by week </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:right;"> 0.028 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 0.028 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:right;"> 0.028 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:right;"> 0.034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 0.034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NFARI probabilities multipliers for X=0 (phi_beta) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NFARI probabilities multipliers for vaccinated or protected persons (theta_beta) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of visit for a case of FARI (delta2): </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> For unvaccinated or unprotected with U=1 </td>
   <td style="text-align:right;"> 0.400 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Multiplier for vaccinated (theta_delta2) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Multiplier for U=0 (mu_delta2) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of visit for a case of NFARI (delta1): </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> For unvaccinated or unprotected with U=1 </td>
   <td style="text-align:right;"> 0.250 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Multiplier for vaccinated (theta_delta1) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Multiplier for U=0 (mu_delta1) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Probability of testing positive for influenza: </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Given influenza infection (tau2) </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Given no influenza infection (tau1) </td>
   <td style="text-align:right;"> 0.000 </td>
  </tr>
</tbody>
</table>

### Running The Code
To run the simulation program using an input CSV file, first run `model4.trueve()` and then using the output from the true VE program, run `model4.flusim()` using the same input file.

```r
  trueve<-model4.trueve("inputfile.csv",vaccov=.495)
  model4.flusim("inputfile.csv",ve.tsi=trueve[3,],ve.tmai=trueve[4,])
```

