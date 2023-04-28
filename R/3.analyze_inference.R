#######################################################################
##                        Analyze MCMC chains 
#######################################################################

rm(list = ls())
library(RColorBrewer)
library(tidyverse)
library(ggpubr)
library(rlang)
library(coda)
library(ggridges)

# Color scheme
cols2 = c("correct" = "royalblue4", "incorrect" = "lightsalmon")

# Model parameters
alpha = 1e-3

# Load post-processing files-------------------------
cond = "flu_covid" 
estimates = read.table(paste0("results/simulation_estimates_", cond, ".txt"), sep = "\t",
                       header = T)
ess = read.table(paste0("results/simulation_ess_", cond, ".txt"), sep = "\t", header = T)
p_accept = read.table(paste0("results/simulation_paccept_", cond, ".txt"), sep = "\t", header = T)
correlations = read.table(paste0("results/simulation_correlations_", cond, ".txt"), sep = "\t", header = T)

# ESS values------------------------------------------
ess %>%
  pivot_longer(c(alpha, beta, rSC, rInfSC), names_to = "param", values_to = "value") %>%
  mutate(param = recode(param, 
                        "alpha" = "community\nfoi",
                        "beta" = "baseline\nhousehold\nfoi", 
                        "rSC" = "child\nrelative\nsusceptibility", 
                        "rInfSC" = "child\nrelative\ninfectivity")) %>%
  ggplot(., aes(x = param, y = value, col = correct_inference)) +
  facet_grid(cols = vars(combination)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3), alpha = 0.3) +
  geom_hline(yintercept = 500) +
  theme_light() +
  labs(col = "Inference", x = "", y = "Effective sample size") +
  scale_color_manual(values = cols2) +
  ylim(c(0,3500))


# Acceptance probability-------------------------------
p_accept %>%
  filter(param != "data") %>%
  mutate(param = recode(param, 
                        "alpha" = "community\nfoi",
                        "beta" = "baseline\nhousehold\nfoi", 
                        "rSC" = "child\nrelative\nsusceptibility", 
                        "rInfSC" = "child\nrelative\ninfectivity")) %>%
  ggplot(., aes(x = param, y = prob, col = correct_inference)) +
  facet_grid(cols = vars(combination)) +
  geom_boxplot() +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3), alpha = 0.3) +
  theme_light() +
  labs(col = "Inference", x = "", y = "Acceptance probability") +
  scale_color_manual(values = cols2) +
  ylim(c(0,1))


# Correlations-------------------------------
correlations %>%
  pivot_longer(c(beta_sus, beta_inf, sus_inf), names_to = "correlation", values_to = "value") %>%
  mutate(correlation = recode(correlation, 
                              "beta_inf" = "foi\ninfectivity",
                              "beta_sus" = "foi\nsusceptibility", 
                              "sus_inf" = "infectivity\nsusceptibility")) %>%
  ggplot(., aes(x = correlation, y = value, col = correct_inference)) +
  facet_grid(cols = vars(combination)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.3), alpha = 0.3) +
  theme_light() +
  labs(col = "Inference", x = "", y = "Correlation") +
  scale_color_manual(values = cols2) +
  ylim(c(-1,1))


# Comparison prior and posterior--------------
flu_pp = list()
post_dist = bind_rows(
  read.table("results/flu/hetero_hetero/mcmc_1.txt", header = T, sep = " ") %>%
    slice(round(0.1*n()):n()) %>%
    mutate(inference = "correct"),
  read.table("results/flu/hetero_homo/mcmc_1.txt", header = T, sep = " ") %>%
    slice(round(0.1*n()):n()) %>%
    mutate(inference = "incorrect")
)

flu_pp[[1]] = post_dist %>%
  ggplot(., aes(x = alpha)) +
  stat_function(aes(col = "prior"), fun = function(x) {dexp(x, 500)}) +
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(alpha), col = "") +
  xlim(c(0,0.01))

flu_pp[[2]] = post_dist %>%
  ggplot(., aes(x = beta)) +
  stat_function(aes(col = "prior"), fun = function(x) {dunif(x, 0, 10)})+
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(frac(beta,4/2)*kappa[child-child]*mu[sus-child]*mu[inf-child]), col = "") +
  xlim(c(0,1))

flu_pp[[3]] = post_dist %>%
  ggplot(., aes(x = rSC)) +
  stat_function(aes(col = "prior"), fun = function(x) {dlnorm(x, 0, 1)})+
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(mu[sus-adult]), col = "") +
  xlim(c(0,3.5))

flu_pp[[4]] = post_dist %>%
  ggplot(., aes(x = rInfSC)) +
  stat_function(aes(col = "prior"), fun = function(x) {dlnorm(x, 0, 1)})+
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(mu[inf-adult]), col = "") +
  xlim(c(0,3.5))

flu_to_save = ggarrange(plotlist = flu_pp, ncol = 4,
                        common.legend = T, legend = "bottom", align = "hv")
ggsave("figures/flu_comparison_prior_posterior.png", flu_to_save, height = 4, width = 11)
ggsave("Figures/flu_comparison_prior_posterior.pdf", flu_to_save, height = 4, width = 11)



covid_pp = list()
post_dist = bind_rows(
  read.table("results/covid/hetero_hetero/mcmc_1.txt", header = T, sep = " ") %>%
    slice(round(0.1*n()):n()) %>%
    mutate(inference = "correct"),
  read.table("results/covid/hetero_homo/mcmc_1.txt", header = T, sep = " ") %>%
    slice(round(0.1*n()):n()) %>%
    mutate(inference = "incorrect")
)

covid_pp[[1]] = post_dist %>%
  ggplot(., aes(x = alpha)) +
  stat_function(aes(col = "prior"), fun = function(x) {dexp(x, 500)}) +
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(alpha), col = "") +
  xlim(c(0,0.01))

covid_pp[[2]] = post_dist %>%
  ggplot(., aes(x = beta)) +
  stat_function(aes(col = "prior"), fun = function(x) {dunif(x, 0, 10)})+
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(frac(beta,4/2)*kappa[child-child]*mu[sus-child]*mu[inf-child]), col = "") +
  xlim(c(0,1))

covid_pp[[3]] = post_dist %>%
  ggplot(., aes(x = rSC)) +
  stat_function(aes(col = "prior"), fun = function(x) {dlnorm(x, 0, 1)})+
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(mu[sus-adult]), col = "") +
  xlim(c(0,3.5))

covid_pp[[4]] = post_dist %>%
  ggplot(., aes(x = rInfSC)) +
  stat_function(aes(col = "prior"), fun = function(x) {dlnorm(x, 0, 1)})+
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(inference~.) +
  scale_color_manual(values = c("posterior" = "orchid4", "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", title = expression(mu[inf-adult]), col = "") +
  xlim(c(0,3.5))

covid_to_save = ggarrange(plotlist = covid_pp, ncol = 4,
                          common.legend = T, legend = "bottom", align = "hv")
ggsave("figures/covid_comparison_prior_posterior.png", covid_to_save, height = 4, width = 11)
ggsave("Figures/covid_comparison_prior_posterior.pdf", covid_to_save, height = 4, width = 11)


# Plot alpha------------------------------------------
a_estimate = estimates %>%
  filter(param == "alpha") %>%
  ggplot(., aes(x = correct_inference, y = median)) +
  facet_grid(rows = vars(combination), scales = "free_y") +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = alpha)) +
  geom_jitter(aes(color = correct_inference)) + 
  theme_light() +
  theme(legend.position = "none") +
  labs(x="", y = "", title = "Alpha")
ggsave(paste0("figures/alpha_",cond,"_estimates.png"), a_estimate, width = 8, height = 6)

alpha_stats = estimates %>%
  filter(param == "alpha") %>%
  group_by(correct_inference, combination) %>%
  summarise(calibration = sum(q2_5<=alpha & alpha<=q97_5)/10,
            error = mean((median-alpha)/alpha)*100,
            .groups = "drop_last")

a_err = alpha_stats %>%
  ggplot(., aes(x = correct_inference, y = error, fill = correct_inference)) +
  facet_grid(col = vars(combination)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values =cols2) +
  theme_light() +
  theme(legend.position = "none") + 
  labs(x = "", y = "Relative bias", col = "", title = "Alpha") 

a_cal = alpha_stats %>%
  ggplot(., aes(x = correct_inference, y = calibration, fill = correct_inference)) +
  facet_grid(col = vars(combination)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values =cols2) +
  theme_light() +
  theme(legend.position = "none") + 
  ylim(c(0,100)) +
  labs(x = "", y = "Calibration", title = "Alpha") 

a = ggarrange(a_err, a_cal, nrow = 2)
ggsave(paste0("figures/alpha_",cond,"_stats.png"), a, width = 5, height = 8)

# Plot foi, rSC, rInfSC--------------------------------------
for (disease in c("flu", "covid")) {
  
  plots = vector("list", 9)
  
  # FOI
  plots[[1]] = estimates %>%
    filter(combination == disease, param == "beta") %>%
    mutate(beta = ifelse(combination == "flu", 2*0.35*1.31/(2*1), 2*0.14*1.31/(0.8*0.5))) %>%
    ggplot(., aes(x = correct_inference, y = median, col = correct_inference)) +
    geom_hline(aes(yintercept = beta)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.3, alpha = 0.3) +
    scale_color_manual(values = cols2) +
    theme_light() +
    theme(legend.position = "none") +
    labs(x = "", y = "Median estimate")
  
  b_df = estimates %>%
    filter(combination == disease, param == "beta") %>%
    mutate(b = ifelse(combination == "flu", 2*0.35*1.31/(2*1), 2*0.14*1.31/(0.8*0.5))) %>%
    group_by(correct_inference, data) %>%
    summarise(
      calibration = sum(q2_5<=b & b<=q97_5)/10,
      error = mean((median-b)/b)*100
    )
  
  plots[[4]] = b_df %>%
    ggplot(., aes(x = correct_inference, y = error, fill = correct_inference)) +
    geom_bar(stat ="identity", width = 0.3) +
    scale_fill_manual(values = cols2) +
    theme_light() +
    theme(legend.position = "none") + 
    ylim(c(-30,30)) +
    labs(x = "", y = "Mean relative bias") 
  
  plots[[7]] = b_df %>%
    ggplot(., aes(x = correct_inference, y = calibration, fill = correct_inference)) +
    geom_bar(stat ="identity", width = 0.3) +
    scale_fill_manual(values = cols2) +
    theme_light() +
    ylim(c(0,100)) +
    theme(legend.position = "none") + 
    labs(x = "", y = "Calibration") 
  
  # rSC and rInfSC 
  i=2
  for (p in c("rSC", "rInfSC")) {
    
    plots[[i]] = estimates %>%
      filter(combination == disease, param == p) %>%
      mutate(rSC = ifelse(combination == "flu", 2, 0.5), 
             rInfSC = ifelse(combination == "flu", 1, 0.8)) %>%
      ggplot(., aes(x = correct_inference, y = median, col = correct_inference)) +
      geom_hline(aes(yintercept = !!sym(p))) +
      geom_boxplot(outlier.shape = NA) +
      scale_color_manual(values=cols2) +
      geom_jitter(width = 0.3, alpha = 0.3) + 
      theme_light() +
      theme(legend.position = "none") +
      labs(x = "", y = "Median estimate")
    
    
    r_df = estimates %>%
      filter(combination == disease, param == p) %>%
      mutate(rSC = ifelse(combination == "flu", 2, 0.5), 
             rInfSC = ifelse(combination == "flu", 1, 0.8)) %>%
      group_by(correct_inference, data) %>%
      summarise(
        calibration = sum(q2_5<=!!sym(p) & !!sym(p)<=q97_5)/10,
        error = mean((median-!!sym(p))/!!sym(p))*100
      )
    
    plots[[i+3]] = r_df %>%
      ggplot(., aes(x = correct_inference, y = error, fill = correct_inference)) +
      geom_bar(stat = "identity", width = 0.3) +
      scale_fill_manual(values = cols2) +
      theme_light() +
      theme(legend.position = "none") + 
      ylim(c(-30,30)) +
      labs(x = "", y = "Mean relative bias") 
    
    plots[[i+6]] = r_df %>%
      ggplot(., aes(x = correct_inference, y = calibration, fill = correct_inference)) +
      geom_bar(stat = "identity", width = 0.3) +
      scale_fill_manual(values = cols2) +
      theme_light() +
      theme(legend.position = "none") + 
      ylim(c(0,100)) +
      labs(x = "", y = "Calibration") 
    
    i=3
  }
  
  tosave = ggarrange(plotlist = plots, ncol = 3, nrow = 3, 
                     labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), align = "hv")
  ggsave(paste0("figures/", disease, ".png"), 
         tosave, height = 6, width = 6)
  ggsave(paste0("figures/", disease, ".pdf"), 
         tosave, height = 6, width = 6)
}


# Summary statistics--------------------------------------
estimates %>%
  filter(param == "beta") %>%
  mutate(p = ifelse(combination == "flu", 2*0.35*1.31/(2*1), 2*0.14*1.31/(0.8*0.5))) %>%
  group_by(correct_inference, combination) %>%
  summarise(
    calibration = sum(q2_5<=p & p<=q97_5)/10,
    error = mean((median-p)/p)*100
  )

estimates %>%
  filter(param == "rSC") %>%
  mutate(rSC = ifelse(combination == "flu", 2, 0.5), 
         rInfSC = ifelse(combination == "flu", 1, 0.8)) %>%
  group_by(correct_inference, param, combination) %>%
  summarise(
    calibration = sum(q2_5<=rSC & rSC<=q97_5)/10,
    error = mean((median-rSC)/rSC)*100
  )

estimates %>%
  filter(param == "rInfSC") %>%
  mutate(rSC = ifelse(combination == "flu", 2, 0.5), 
         rInfSC = ifelse(combination == "flu", 1, 0.8)) %>%
  group_by(correct_inference, param, combination) %>%
  summarise(
    calibration = sum(q2_5<=rInfSC & rInfSC<=q97_5)/10,
    error = mean((median-rInfSC)/rInfSC)*100
  )

estimates %>%
  filter(param == "pInf") %>%
  mutate(p = ifelse(combination == "flu", 1-exp(-2*0.35*1.31/(2*1)), 1-exp(-2*0.14*1.31/(0.8*0.5)))) %>%
  group_by(correct_inference, combination) %>%
  summarise(
    calibration = sum(q2_5<=p & p<=q97_5)/10,
    error = mean((median-p)/p)*100
  )

