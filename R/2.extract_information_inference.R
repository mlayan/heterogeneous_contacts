#######################################################################
##                        Analyze MCMC chains 
#######################################################################

rm(list = ls())
library(tidyverse)
library(coda)

# Parameter combination 
param_comb = c("flu", "covid")

# Post-processing of posterior distributions
estimates = data.frame()
p_accept = data.frame()
ess = data.frame()
correlations = data.frame()

for (condition in c("hetero_hetero", "hetero_homo")) {
  for (comb in param_comb) {
    for (f in list.files(paste0("results/", comb, "/", condition), pattern = "mcmc_[0-9]*.txt$",
                         full.names = T)) {

      I = paste0(gsub(".*_", "", condition), "geneous")
      D = paste0(gsub("_.*", "", condition), "geneous")
      post_tmp = read.table(f, header = T, sep=" ")
      post_tmp$nSim = as.numeric(gsub(".*_|.txt", "", f))
      post_tmp$inference = I
      post_tmp$data = D
      post_tmp$combination = comb
      post_tmp$correct_inference = ifelse(I == D, "correct", "incorrect")
      post_tmp = post_tmp[ceiling(0.1*nrow(post_tmp)):nrow(post_tmp),]
      post_tmp$rSC = 1/post_tmp$rSC
      post_tmp$rInfSC = 1/post_tmp$rInfSC
      post_tmp$beta = 2*post_tmp$beta/(post_tmp$rSC * post_tmp$rInfSC * 0.76)
      post_tmp$pInf = 1-exp(-post_tmp$beta)

      # Estimates of the posterior distribution
      estimates_tmp = post_tmp %>%
        select(nSim, inference, data, correct_inference, combination, alpha, beta, rSC, rInfSC, pInf) %>%
        pivot_longer(cols = c(alpha, beta, rSC, rInfSC, pInf), values_to = "value", names_to = "param") %>%
        group_by(nSim, inference, data, correct_inference, combination, param) %>%
        summarise(median = median(value),
                  q97_5 = quantile(value, 0.975),
                  q2_5 = quantile(value,0.025), .groups = "drop_last")
      estimates = bind_rows(estimates, estimates_tmp)

      # Probability of acceptance
      p_accept_tmp = post_tmp %>%
        select(., matches("_[a,p]$")) %>%
        summarise(across(where(is.numeric), ~ sum(.x))) %>%
        pivot_longer(everything()) %>%
        mutate(param = gsub("_.*$", "", name), type = gsub("^.*_", "", name)) %>%
        select(param, type, value) %>%
        pivot_wider(values_from = value, names_from = type) %>%
        filter(p > 0) %>%
        mutate(prob = a/p) %>%
        select(param, prob) %>%
        mutate(nSim = as.numeric(gsub(".*_|.txt", "", f)), combination = comb,
               correct_inference = ifelse(I == D, "correct", "incorrect"))
      p_accept = bind_rows(p_accept, p_accept_tmp)

      # Chain mixing
      ess_tmp = post_tmp %>%
        select(alpha, beta, rSC, rInfSC) %>%
        summarise(across(where(is.numeric), ~effectiveSize(.x))) %>%
        mutate(nSim = as.numeric(gsub(".*_|.txt", "", f)), combination = comb,
               correct_inference = ifelse(I == D, "correct", "incorrect"))
      ess = bind_rows(ess, ess_tmp)

      # Correlations
      correlations_tmp = data.frame(
        beta_sus=cor(post_tmp$beta, post_tmp$rSC),
        beta_inf=cor(post_tmp$beta, post_tmp$rInfSC),
        sus_inf=cor(post_tmp$rSC, post_tmp$rInfSC),
        nSim=as.numeric(gsub(".*_|.txt", "", f)),
        correct_inference=ifelse(I == D, "correct", "incorrect"),
        combination=comb
        )
      correlations = bind_rows(correlations, correlations_tmp)

    }
  }
}

write.table(estimates, "results/simulation_estimates_flu_covid.txt", sep = "\t",
            quote = F, row.names = F)

write.table(ess, "results/simulation_ess_flu_covid.txt", sep = "\t",
            quote = F, row.names = F)

write.table(p_accept, "results/simulation_paccept_flu_covid.txt", sep = "\t",
            quote = F, row.names = F)

write.table(correlations, "results/simulation_correlations_flu_covid.txt", sep = "\t",
            quote = F, row.names = F)


