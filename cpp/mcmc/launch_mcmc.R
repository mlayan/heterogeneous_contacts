#!/usr/bin/env Rscript

# Script header
header <- "#!/bin/sh
#SBATCH --mem=500
#SBATCH -p mmmi -q mmmi"

# Script footer
footer <- "exit 0"

# Parameters
library(dplyr, quietly = T)
#param_comb = apply(expand.grid(seq(0.5, 2, 0.25), # rSC
#                    seq(0.5, 2, 0.25), # rInfSC
#                    stringsAsFactors = F),
#                1, paste, collapse = "_")
#param_comb = c("1.75_1")

param_comb = c("flu", "covid")

params = expand.grid(1:1000, 
                     "heterogeneous", # data
                     c("heterogeneous","homogeneous"), # inference
                     param_comb,
                     stringsAsFactors = F)

for (curr in 1:nrow(params)) {
  # Redirect output and error messages
  reOut <- paste0("#SBATCH -o mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".log")
  reErr <- paste0("#SBATCH -e mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".err")
  
  # Command to run the script
  command <- paste("srun program.out", paste0(t(params[curr, ]), collapse = " "))
  
  # Create script to submit
  scriptToLaunch <- paste0("mcmc_", paste0(t(params[curr, ]), collapse = "_"), ".sh")
  
  script <- paste0(header, "\n",
                   reOut, "\n", 
                   reErr, "\n\n",
                   command, "\n\n",
                   paste0("rm ", scriptToLaunch, "\n"),
                   footer)
  
  
  
  write(script, scriptToLaunch)
  
  # Submit to cluster queue
  #system(paste("sbatch", scriptToLaunch))
}
