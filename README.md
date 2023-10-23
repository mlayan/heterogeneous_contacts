# Addressing current limitations of household transmission studies by collecting contact data

Full text is available [here]().

## Description

Here, we investigate the impact of contact heterogeneity within households on the estimation of the relative infectivity and relative susceptibility of children compared to adults.  

### Simulating household transmission chains  

In `R/1.simulate_epidemics.R`, one can simulate household epidemics assuming heterogeneous contact patterns using a discretized transmission model coded in Rcpp and available in `cpp/simulation/simulate_epidemic.cpp`. The *in silico* household data set is derived from the RECOVER study available [here](https://doi.org/10.1007/s10654-022-00870-9). Of note, mothers are coded as 0, fathers as 1, and children as 2.

If you want to modify the length of the study period, you can specify it in the R code. 

Examples of 20 simulation outputs are available in `data` for both the flu and the COVID-19 scenarios.

### Inference

In `cpp/mcmc/`, one can find the inference model. The R script `launch_mcmc.R` allows to launch the MCMC chains for all simulated epidemics with the homogeneous or heterogeneous mixing assumption.

To compile the code, you need a gcc compiler (version 9.3.0 or 9.2.0) and the boost library (version 1.72.0).

The current version of the model estimates the force of infection in households of size 4. 

Examples of inference outputs are available in `results` for both the flu and the COVID-19 scenarios.

### Analysis of inference chains

To analyze the inferred parameters for each simulated epidemic, household mixing assumption (homogeneous or heterogeneous), and disease scenario, one can first extract the summary statistics of MCMC chains using the `R/2.extract_information_inference.R` code. This scripts loads all the MCMC chains and computes effective sampling sizes (ESS), probability of acceptance, summary statistics of the post burn-in chains, and correlations between parameters. These information are stored in the `results` directory. 

Finally, figures are produced using `R/3.analyze_inference.R` and stored in `figures`.

