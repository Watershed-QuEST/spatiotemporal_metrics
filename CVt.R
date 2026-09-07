#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-08-23 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-08

# Requires: 02_build_synthetic_data.R must be run first to produce data/nm_clean.csv, data/nm_field_setups.rds, and data/nm_synthetic_extended.csv. If these are available, no need to rerun 02_build_synthetic_data.R

# This script demonstrates how to calculate CVt, the temporal coefficient of variation (across sampling campaigns, one value per sampling site). 
# It also performs a sensitivity analysis to answer the question: how much does a CVt estimate vary/undershoot depending how many sampling campaigns are done? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possible fits to a SSN2 spatial covariance model. 
# Ref_CVt (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R, which simulates 60 (1 per month x 5 years) synoptic campaigns.

# This script is organized as:
#   PART A -- CVt_observed: computed directly from the real nm data
#   PART B -- Monte Carlo sensitivity analysis

# Outputs:
#   data/CVt_observed.csv:      real CVt per site x constituent
#   data/CVt_observed_mean.csv: real CVt averaged across sites, per constituent
#   data/CVt_MC_summary.csv:    sensitivity analysis summary (median, P05-P95, bias vs Ref_CVt)
#   plots/CVt_MC_plot.png

#### Packages ####
library(tidyverse)
#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

# Monte Carlo iterations
n_iter    <- 10000
camp_grid <- c(2, 3, 4, 5, 7, 10, 15, 20, 30, 50)      # campaign counts to test

n_synthetic_months <- 60      # must match 02_build_synthetic_data.R
time_unit_days     <- 30.44   # must match 02_build_synthetic_data.R
total_days <- n_synthetic_months * time_unit_days  # span campaign times are drawn from

set.seed(42)

clean              <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups       <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)

actual_n_campaigns <- n_distinct(clean$CampaignID)  # real number of nm campaigns conducted


#### PART A -- Calculate CVt of real toy dataset ####
# Per-site CV across that site's real campaigns

## For each constituent in `constituents`: per-site CV across that site's real campaigns (sites with fewer than 3 real campaigns are dropped, since sd() of 1-2 points is too noisy to call a temporal CV).
CVt_observed_list      <- list()
CVt_observed_mean_list <- list()

for (cc in constituents) {
  
  per_site_cvt <- clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(Site) %>%
    summarize(n_campaigns = n(), CVt = sd(.data[[cc]]) / mean(.data[[cc]]), .groups = "drop") %>%
    filter(n_campaigns >= 3)

  print(per_site_cvt)
  
  cvt_mean <- mean(per_site_cvt$CVt)
  
  CVt_observed_list[[cc]]      <- per_site_cvt %>% mutate(Constituent = cc, .before = 1)
  CVt_observed_mean_list[[cc]] <- tibble(Constituent = cc, n_sites_used = nrow(per_site_cvt),
                                         CVt_observed = cvt_mean)
}

## ---- Combine and save ----
CVt_observed <- bind_rows(CVt_observed_list)
write_csv(CVt_observed, file.path(data_out_dir, "CVt_observed.csv"))

CVt_observed_mean <- bind_rows(CVt_observed_mean_list)
write_csv(CVt_observed_mean, file.path(data_out_dir, "CVt_observed_mean.csv"))
print(CVt_observed_mean)


#### PART B -- Monte Carlo temporal sensitivity analysis ####

# For each constituent: draw n_iter independent synthetic draws at each campaign count in camp_grid, and see how the resulting CVt estimate's median and spread compare to Ref_CVt (the "full 5-year record" truth from the static extended dataset).

# Each draw does the same four things:
# 1. pick n random campaign dates across the 5-year synthetic span
# 2. build spatially-correlated random noise for those n campaigns (using SSN2's fit, L_space below)
# 3.correlate that noise across the n campaigns (red noise), based on how many days apart each pair of them is, using the fitted persistence estimate (rho)
# 4. add SSN2's spatial mean and seasonal cycle, exponentiate back to the natural (non-log) scale, and compute the per-site CV across those n campaigns

# The loops below iterate over constituents, campaignn, and n_iter random draws 

temporal_mc_results_list <- list()

for (cc in constituents) {
  
  fs       <- field_setups[[cc]]
  m        <- fs$setup$m           # number of predpts locations
  mu_cond  <- fs$setup$mu_cond     # SSN2's kriged mean, one value per location
  L_space  <- fs$setup$L_space     # SSN2's spatial Cholesky factor
  beta     <- fs$beta              # seasonal-cycle coefficients (Intercept, logArea, sin, cos)
  ref_doy  <- fs$ref_doy
  rho      <- fs$temporal_rho      # red-noise persistence, from 02's residual diagnostic
  
  if (rho < 0 || rho >= 1) {
    stop(cc, " temporal_rho (", rho, ") is outside [0, 1) -- this exponential-decay ",
         "temporal correlation can't represent a negative value; see 02_build_synthetic_data.R's Part D.")
  }
  
  # Reference CVt: the static extended dataset's per-site CV across all 60 synthetic months, averaged across all m predpts locations. This is a large, stable estimate of the "full temporal record" "truth"
  ref_cvt <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(predID) %>%
    summarize(CVt = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVt) %>% mean()
  
  # The seasonal shift, evaluated at ref_doy 
  ref_shift <- beta[3] * sin(2 * pi * ref_doy / 365) + beta[4] * cos(2 * pi * ref_doy / 365)
  
  mc_results <- list()
  
  for (n in camp_grid) {
    
    cvt_draws <- numeric(n_iter)
    
    for (i in seq_len(n_iter)) {
      
      # Step 1: draw n random campaign dates across the 5-year synthetic span
      t_sample   <- sort(runif(n, 0, total_days))
      doy_sample <- t_sample %% 365.25
      
      # Step 2: generate spatially-correlated noise, one column per campaign
      z             <- matrix(rnorm(m * n), nrow = m, ncol = n)
      spatial_noise <- t(L_space) %*% z
      
      # Step 3: also correlate those n columns with each other, based on elapsed time between campaigns (red noise, if rho > 0; this reduces to no correlation at all if rho == 0).
      elapsed_days  <- abs(outer(t_sample, t_sample, "-")) / time_unit_days
      time_corr     <- rho ^ elapsed_days
      time_chol     <- chol(time_corr + diag(1e-8, n))
      noise         <- spatial_noise %*% t(time_chol)
      
      # Step 4: deterministic seasonal shift for these campaign dates (relative to ref_doy), then the full synthetic field: SSN2's spatial mean + seasonal shift + noise, exponentiated back to the natural (non-log) scale.
      shift     <- beta[3] * sin(2 * pi * doy_sample / 365) +
        beta[4] * cos(2 * pi * doy_sample / 365) - ref_shift
      shift_mat <- matrix(shift, nrow = m, ncol = n, byrow = TRUE)
      field     <- exp(mu_cond + shift_mat + noise)  # m x n, natural scale
      
      # Per-site CV across these n campaigns, averaged across all locations. 
      row_mean     <- rowMeans(field)
      row_sd       <- sqrt(rowSums((field - row_mean)^2) / (n - 1))
      site_cvt     <- row_sd / row_mean
      cvt_draws[i] <- mean(site_cvt)
    }
    
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n, CVt = cvt_draws, Ref_CVt = ref_cvt)
  }
  
  temporal_mc_results_list[[cc]] <- bind_rows(mc_results)
}

## Combine all constituents, summarize, and plot

temporal_mc_results <- bind_rows(temporal_mc_results_list)

temporal_mc_summary <- temporal_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVt = median(CVt), P05_CVt = quantile(CVt, 0.05), P95_CVt = quantile(CVt, 0.95),
            Ref_CVt = first(Ref_CVt), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVt - Ref_CVt) / Ref_CVt * 100)

print(temporal_mc_summary, n = Inf)

write_csv(temporal_mc_summary, file.path(data_out_dir, "CVt_MC_summary.csv"))

# One row per constituent, the analysis's own Pct_Bias at N = actual_n_campaigns
bias_annot <- temporal_mc_summary %>%
  filter(N == actual_n_campaigns) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias", N, Pct_Bias))

p_cvt <- ggplot(temporal_mc_summary, aes(x = N, y = Median_CVt)) +
  geom_ribbon(aes(ymin = P05_CVt, ymax = P95_CVt), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVt), color = "blue", linetype = "dashed") +
  geom_vline(xintercept = actual_n_campaigns, color = "red", linetype = "dashed") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N campaigns)", y = "CVt (median, 5th-95th pct band)",
       title = "Temporal sample-size sensitivity",
       subtitle = "Blue = reference CVt (static extended dataset, all 60 synthetic months); red = actual sample size, annotated with pct bias at that N. Each Monte Carlo iteration redraws a fresh spatiotemporal field, campaign times drawn at random across a 5-year span.") +
  theme_bw()
ggsave(file.path(plot_dir, "CVt_MC_plot.png"), p_cvt, width = 8, height = 6, dpi = 150)

