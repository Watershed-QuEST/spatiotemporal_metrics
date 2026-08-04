#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-08-03

# Requires: 02_build_synthetic_data.R must be run first (produces data/nm_clean.csv, data/nm_field_setups.rds, and data/nm_synthetic_extended.csv), and spatiotemporal_helpers.R must be in the same folder as this script.

# This script demonstrates how to calculate CVs, the spatial coefficient of variation (across sites, one value per sampling event). 
# It also performs a sensitivity analysis to answer the question: how much does a CVs estimate vary/undershoot depending on how many sites you sample? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possibile fits to a SSN2 spatial covariance model. 
# Ref_CVs (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R, which simulates a synoptic campaign that sampled every 250  m along the network.

# This script is organized as:
#   PART A -- CVs_observed: computed directly from the real nm data
#   PART B -- Monte Carlo sensitivity sweep: 

# To adapt this script to a metric other than CVs: swap out the statistic computed inside the two n_iter loops in Part B, and Part A's real-data calculation to match.

# Outputs:
#   data/CVs_by_campaign.csv:  real CVs, one row per real campaign x constituent
#   data/CVs_observed.csv: Mean_CVs and SD_CVs across those real campaigns
#   data/CVs_MC_summary.csv:  MC summary
#   plots/CVs_MC_plot.png: CVs vs N, with reference line and bias annotation
#   plots/CVs_sensitivity_metrics_plot.png: percent bias and SD vs N, side by side

#### Packages ####
library(tidyverse)
source("spatiotemporal_helpers.R")

#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

n_iter    <- 10000                                  # Monte Carlo iterations
site_grid <- c(3:50, 75, 100, 150, 300)              # site counts to test

actual_n_sites <- 23  # real number of nm sites sampled -- update if that changes

set.seed(42)

clean            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups     <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)

#### PART A -- Calculate CVs of real toy dataset ####
# Each site's concentration averaged across its available real campaigns, then CV taken across sites.

CVs_by_campaign <- map_dfr(names(field_setups), function(cc) {
  clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(CampaignID) %>%
    summarize(n_sites = n(),
              CVs = if (n() >= 2) sd(.data[[cc]]) / mean(.data[[cc]]) else NA_real_,
              .groups = "drop") %>%
    mutate(Constituent = cc, .before = 1)
})

print(CVs_by_campaign)
write_csv(CVs_by_campaign, file.path(data_out_dir, "CVs_by_campaign.csv"))

## Summarized across campaigns per constituent
CVs_observed <- CVs_by_campaign %>%
  filter(!is.na(CVs)) %>%
  group_by(Constituent) %>%
  summarize(n_campaigns = n(), Mean_CVs = mean(CVs),
            SD_CVs = if (n() >= 2) sd(CVs) else NA_real_, .groups = "drop")

print(CVs_observed)
write_csv(CVs_observed, file.path(data_out_dir, "CVs_observed.csv"))

#### PART B -- Monte Carlo spatial sensitivity analysis ####

spatial_mc_results <- list()

for (cc in names(field_setups)) {
  fs <- field_setups[[cc]]
  m <- fs$setup$m
  this_site_grid <- site_grid[site_grid <= m]
  if (length(this_site_grid) == 0) {
    cat("  SKIP", cc, "-- every site_grid value exceeds m =", m, "predpts.\n")
    next
  }
  
  cat("\n---", cc, "[spatial MC,", n_iter, "iterations, up to", max(this_site_grid), "sites] ---\n")
  
  ## Calc reference CVs
  ref_cvs <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  ## One fresh, spatially-correlated noise draw per Monte Carlo iteration; drawn all at once (an m x n_iter matrix in one call). t_days is just a placeholder vector of the right length.
  noise_pool <- draw_field(fs$setup, t_days = seq_len(n_iter), rho = 0)
  sim_mat <- exp(fs$setup$mu_cond + noise_pool)  # m x n_iter, natural scale
  
  for (n in this_site_grid) {
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx <- sample(m, n)
      vals <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    spatial_mc_results[[paste(cc, n)]] <- tibble(Constituent = cc, N = n, CVs = cvs_draws, Ref_CVs = ref_cvs)
  }
}

spatial_mc_results <- bind_rows(spatial_mc_results)

# present and save results
if (nrow(spatial_mc_results) == 0) {
  cat("\nspatial_mc_results is empty -- see the SKIP message(s) above.\n")
} else {
  spatial_mc_summary <- spatial_mc_results %>%
    group_by(Constituent, N) %>%
    summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
              SD_CVs = sd(CVs), Ref_CVs = first(Ref_CVs), .groups = "drop") %>%
    mutate(Pct_Bias = (Median_CVs - Ref_CVs) / Ref_CVs * 100)
  
  cat("\n--- Spatial (CVs) sample-size sensitivity ---\n")
  print(spatial_mc_summary, n = Inf)
  
  write_csv(spatial_mc_summary, file.path(data_out_dir, "CVs_MC_summary.csv"))
  
  ## Creates one row per constituent
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
  
  ## Second plot: percent bias and precision (SD across Monte Carlo iterations), both as a function of N. Pct_Bias shows how far off the median estimate tends to be; SD_CVs shows how much that estimate itself changes at a given sample size. Pivoted to long format so both metrics can share one facet grid (metric x constituent) rather than needing a dual y-axis.
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
  
}
