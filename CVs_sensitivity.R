#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-09

# Requires: 02_build_synthetic_data.R must be run first (produces data/*_clean.csv, data/*_field_setups.rds, and data/*_synthetic_extended.csv).

# This script performs a sensitivity analysis to answer the question: how much does a CVs estimate vary/undershoot depending on how many sites you sample? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possibile fits to a SSN2 spatial covariance model. 
# Ref_CVs (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R, which simulates a synoptic campaign that sampled every 250 m along the network.

# Outputs:
#   data/CVs_MC_summary.csv:  MC summary
#   plots/CVs_MC_plot.png: CVs vs N, with reference line and bias annotation
#   plots/CVs_sensitivity_metrics_plot.png: percent bias and SD vs N, side by side

#### Packages ####
library(tidyverse)

#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_iter    <- 10000                                  # Monte Carlo iterations
site_grid <- c(3:50, 75, 100, 150, 300)              # site counts to test

actual_n_sites <- 23  # real number of nm sites sampled -- update if that changes

set.seed(42)

clean            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups     <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)


#### Monte Carlo spatial sensitivity analysis ####

spatial_mc_results_list <- list()

for (cc in constituents) {
  
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m           # number of predpts locations
  mu_cond <- fs$setup$mu_cond     # SSN2's kriged mean, one value per location
  L_space <- fs$setup$L_space     # SSN2's spatial Cholesky factor
  
  this_site_grid <- site_grid[site_grid <= m]
  
  ## Reference CVs
  ref_cvs <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  ## Step 1-2: one fresh, spatially-correlated snapshot per Monte Carlo iteration, drawn all at once (an m x n_iter matrix in one call) for  speed, which is equivalent to drawing them one at a time: each column is still an independent draw. Independent random numbers, "smeared" by the spatial Cholesky factor so nearby predpts locations end up correlated the way SSN2 fit them.
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)  # m x n_iter, natural scale
  
  mc_results <- list()
  
  for (n in this_site_grid) {
    
    cvs_draws <- numeric(n_iter)
    
    for (i in seq_len(n_iter)) {
      # Step 3: N of the m locations, from this iteration's own snapshot
      idx          <- sample(m, n)
      vals         <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n, CVs = cvs_draws, Ref_CVs = ref_cvs)
  }
  
  spatial_mc_results_list[[cc]] <- bind_rows(mc_results)
}

## Combine all constituents, summarize, and plot

spatial_mc_results <- bind_rows(spatial_mc_results_list)

spatial_mc_summary <- spatial_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs = first(Ref_CVs), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVs - Ref_CVs) / Ref_CVs * 100)

write_csv(spatial_mc_summary, file.path(data_out_dir, "CVs_MC_summary.csv"))

bias_annot <- spatial_mc_summary %>%
  filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias", N, Pct_Bias))

p_cvs <- ggplot(spatial_mc_summary, aes(x = N, y = Median_CVs)) +
  geom_ribbon(aes(ymin = P05_CVs, ymax = P95_CVs), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVs), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dashed") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites)", y = "CVs (median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity",
       subtitle = "Blue = reference CVs; red = actual sample size") +
  theme_bw()
ggsave(file.path(plot_dir, "CVs_MC_plot.png"), p_cvs, width = 8, height = 6, dpi = 150)

## Second plot: percent bias and precision (SD across Monte Carlo iterations)
sensitivity_long <- spatial_mc_summary %>%
  select(Constituent, N, Pct_Bias, SD_CVs) %>%
  pivot_longer(c(Pct_Bias, SD_CVs), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, Pct_Bias = "Percent bias vs Ref_CVs (%)",
                         SD_CVs = "SD of CVs across MC iterations"))

p_sensitivity <- ggplot(sensitivity_long, aes(x = N, y = Value)) +
  geom_line() + geom_point(size = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dashed") +
  facet_grid(Metric ~ Constituent, scales = "free_y") +
  labs(x = "Sample size (N sites)", y = NULL,
       title = "Spatial sensitivity: bias and precision vs. sample size",
       subtitle = "Red = actual sample size") +
  theme_bw()
ggsave(file.path(plot_dir, "CVs_sensitivity_metrics_plot.png"), p_sensitivity, width = 9, height = 6, dpi = 150)
