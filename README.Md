# Introduction

## Model Description
We present a dynamic model consisting of five steps. Below we define the model steps, the associated variables (Table 2), and the probabilities determining each variable’s distribution (Table 3). All variables are defined for each member of the study population, and we allow some variables to change over time (measured in weeks). Figure 1 illustrates the possible sources of confounding and bias present in studies designed to evaluate influenza VE [15,16]. Model assumptions are shown in Web Table 1.

**Step 1:** Covariates. We assume that people within the population can be classified with a health status (X) of either “healthy” or “frail” and a health awareness (U) of either “high” or “low”.

**Step 2:** Vaccination. We consider the vaccination scenario where an individual is considered vaccinated if they received the vaccine at least 14 days prior to the study onset (V = 1) or remains unvaccinated throughout the study (V = 0).

**Step 3:** Influenza and non-influenza ARI. During the influenza season, a person may become infected with an influenza virus and develop influenza ARI. Regardless of influenza infection, a person may develop one or more non-influenza ARIs. We define a variable Yj for the illness/infection status in week j as follows: Yj = 0 for no ARI, Yj = 1 for non-influenza ARI, and Yj = 2 for influenza ARI. If a person has both non-influenza ARI and influenza ARI in the same week, we consider them as influenza ARI (i.e., Yj = 2). The distribution of Yj may depend on the person’s vaccination (V) and health (X) status.

**Step 4:** Seeking medical care for ARI. A person with an ARI in week j may seek medical care (Mj). The probability of seeking medical care depends on Yj, as only those individuals who have an ARI may seek medical care, and it may be different for influenza ARI and non-influenza ARI patients. This probability may also depend on V and U.

**Step 5:** Testing for influenza infection. We assume that each person who seeks medical care for ARI is tested for influenza infection. Let Tj denote the binary test result, where T_j=1 or 0 for positive or negative, respectively.

# Simulation Program

## Input File

## Running The Code
To run the simulation program using an input CSV file, first run `model4.trueve()` and then using the output from the true VE program, run `model4.flusim()` using the same input file.
```{r eval=FALSE}
  trueve <- model4.trueve("inputfile.csv",vaccov=.495)
  model4.flusim("inputfile.csv",ve.tsi=trueve[3,],ve.tmai=trueve[4,])
```

