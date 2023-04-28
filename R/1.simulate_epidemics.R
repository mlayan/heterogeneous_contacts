#######################################################################
##                    Simulate household epidemics
#######################################################################

rm(list = ls())
library(Rcpp)
library(tidyverse)

# Load household dataset
load("data/bdd.RData")

#######################################################################
# Simulation model
#######################################################################
# Remove compiled files
filesToRemove = list.files(
  path = "cpp/simulation/",
  pattern = "^.*\\.o$",
  full.names = T
)
if (length(filesToRemove)) sapply(filesToRemove, file.remove)

# Source cpp file
sourceCpp("cpp/simulation/simulate_epidemic.cpp")

#######################################################################
# Simulate epidemics
#######################################################################
pInf = c()
contact_pattern = "heterogeneous"
alpha = 1e-3 # force of infection from the community
nSimulations = 1000 # number of simulated epidemics
followUp = 20 # 20 days since symptom onset of index case
disease = "covid" # "flu"

cat(paste0(disease, " mixing:\n"))

if (grepl("flu", disease)) {
  rSC = 2 # relative susceptibility of children
  rInfSC = 1  # relative infectivity of children
  beta4 = 0.35 # force of infection between two children in a household of size 4
}

if (grepl("covid", disease)) {
  rSC = 0.5 # relative susceptibility of children
  rInfSC = 0.8 # relative infectivity of children
  beta4 = 0.14 # force of infection between two children in a household of size 4
}

if (!dir.exists(paste0("data/", disease, "/", contact_pattern))) {
  dir.create(paste0("data/", disease, "/", contact_pattern), recursive = T)
}

for (i in 1:nSimulations) {
  
  out = lapply(
    bdd, 
    hhEpidemic,
    alpha = alpha,
    beta = beta4,
    rSC = rSC, 
    rInfC= rInfSC,
    followUp = followUp,
    mainHHSize = 4.0, 
    contact_pattern = contact_pattern,
    dt = 0.01 # time step for the simulation 
    ) 
  
  do.call("rbind", out) %>%
    mutate(dds = trunc(dds)) %>%
    mutate(
      infectionStatus = ifelse(dds > studyPeriod & dds != 1000 & index == 0, 0, infectionStatus),
      dds = ifelse(dds > studyPeriod & dds != 1000 & index == 0, 1000, dds)
    ) %>%
    select(indid, hhid, hhsize, dds, infectionStatus, startFollowUp, studyPeriod, adult, ddi) %>%
    write.table(
      paste0("data/", disease, "/", contact_pattern, "/data_", i, ".txt"),
      row.names = F, col.names = F, sep = " ", quote = F
    )
  
  out = do.call("rbind", out) %>%
    mutate(
      infectionStatus = ifelse(dds > studyPeriod & dds != 1000 & index == 0, 0, infectionStatus),
      dds = ifelse(dds > studyPeriod & dds != 1000 & index == 0, 1000, dds)
    )
  
  pInf = c(
    pInf,
    sum(out$ddi != 1000 & out$index == 0) / sum(out$index == 0)
  )

  print(i)
  print(sum(out$ddi != 1000 & out$index == 0) / sum(out$index == 0))
} 

if (disease == "flu") cat("\n\nFlu\n")
if (disease == "flu") cat("\n\nCOVID\n")
summary(pInf)
