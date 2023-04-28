# Contact heterogeneity in households

## Description

Here, we investigate the impact of contact heterogeneity within households on the estimation of the relative infectivity and relative susceptibility of children compared to adults.  

### Simulating household transmission chains  

In `R/5.inference_troubleshooting.R`, one can simulate household epidemics using a discretized transmission model coded in Rcpp and available in `cpp/simulation/simulate_epidemic_performance.cpp`. Of note, symptomatic individuals are coded as 1 and asymptomatic individuals as 2. Children under 12 are coded as 0, teenagers between 12 and 18 are coded as 1, and adults above 18 are coded as 2. 

If you want to modify the length of the study period, you should specify it both in the R code and in the cpp code. 

Examples of simulation outputs are available in `data`.

### Inference

In `cpp/mcmc/`, one can find the inference model. Importantly, individuals with a negative infection status are considered immune and do not contribute to the likelihood (not applicable in the simulation analysis).

To compile the code, you need a gcc compiler (9.3.0 or 9.2.0) and the boost library (1.72.0).

The current version of the model estimates the force of infection in households of size 5 using the true continuous infection and symptom onset dates with no augmentation data step. 

Examples of inference outputs are available in `results`.

One can analyze the output chains using the second part of the `R/5.inference_troubleshooting.R` code.