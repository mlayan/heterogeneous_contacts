#!/usr/bin/env Rscript

# Script header
header <- "#!/bin/sh"

# Script footer
footer <- "exit 0"

# Parameters
library(dplyr, quietly = T)
param_comb = c("flu", "covid")

params = expand.grid(1:1000, 
                     "heterogeneous", # data
                     c("heterogeneous","homogeneous"), # inference
                     param_comb,
                     stringsAsFactors = F)

for (curr in 1:nrow(params)) {
  # Command to run the script
  command <- paste("program.out", paste0(t(params[curr, ]), collapse = " "))
  
  # Create script to submit
  scriptToLaunch <- paste0("mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".sh")
  
  script <- paste0(header, "\n",
                   command, "\n\n",
                   footer)
  
  write(script, scriptToLaunch)
  
  # Submit to cluster queue
  system(paste("sbatch", scriptToLaunch))
}
