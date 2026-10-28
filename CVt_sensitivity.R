#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-10-28

# Requires: 02_build_synthetic_data.R must be run first (produces data/*_clean.csv and data/*_synthetic_extended.csv). Does NOT use *_field_setups.rds -- unlike the spatial version, there's no fitted temporal covariance model to draw fresh realizations from, so this reuses the already-expanded synthetic dataset (all ~60 synthetic campaigns per location) directly as the time-sampling population.

# This script performs a sensitivity analysis to answer the question: how much does a CVt estimate vary/undershoot depending on how many campaigns you sample? Each Monte Carlo iteration picks one synthetic-network location at random (that location's full synthetic time series stands in for "one hypothetical realization," playing the role a fresh covariance-model draw played in the spatial version) and subsamples N campaigns from that location's time series.
# Ref_CVt_MC: the full-time-population CVt, averaged across all synthetic-network locations (internally consistent with the sampling procedure). Ref_CVt_Real: CVt at just the real-equivalent site locations, across all 60 synthetic campaigns.

# Outputs:
#   data/CVt_MC_summary_nm.csv / _br.csv
#   plots/CVt_MC_plot_nm.png / _br.png

#### Packages ####
library(tidyverse)
library(patchwork)

#### ----- nm ----- ####
#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_iter        <- 10000               # Monte Carlo iterations
campaign_grid <- c(3:30, 40, 50, 60)  # campaign counts to test

actual_n_campaigns <- 3  # real number of campaigns sampled -- same for nm and br

set.seed(42)

clean                   <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
synthetic_extended_full <- read_csv(file.path(data_out_dir, "nm_synthetic_with_real.csv"), show_col_types = FALSE)

real_predids <- synthetic_extended_full %>%
  filter(Source == "real") %>%
  pull(predID) %>%
  unique()

synthetic_extended_real_locs <- synthetic_extended_full %>%
  filter(predID %in% real_predids)

#### Varying sample size in time ####

temporal_mc_results_list <- list()
for (cc in constituents) {
  ## Population matrix: rows = synthetic campaigns, columns = locations.
  ## Each column is one location's full realized time series -- the unit we
  ## resample a random one of, per iteration, in place of a fresh MC draw.
  sim_mat_time <- synthetic_extended_full %>%
    filter(Constituent == cc) %>%
    select(CampaignNum, predID, Conc) %>%
    pivot_wider(names_from = predID, values_from = Conc) %>%
    arrange(CampaignNum) %>%
    select(-CampaignNum) %>%
    as.matrix()
  
  m           <- nrow(sim_mat_time)  # number of campaigns in the full synthetic time series
  n_locations <- ncol(sim_mat_time)  # number of locations available to draw from
  this_campaign_grid <- sort(unique(c(campaign_grid[campaign_grid <= m], m)))
  
  ## "Realized truth": CVt at the real-equivalent site locations, computed
  ## over all m campaigns. Tied to the specific locations actually sampled.
  ref_cvt_real <- synthetic_extended_real_locs %>%
    filter(Constituent == cc) %>%
    group_by(predID) %>%
    summarize(CVt = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVt) %>% mean()
  
  ## "Model truth": the full-time-population CVt (all m campaigns), averaged
  ## across every location. Internally consistent with the subsampling loop
  ## below, since it's drawn from the same sim_mat_time.
  full_pop_cvt <- apply(sim_mat_time, 2, function(vals) sd(vals) / mean(vals))
  ref_cvt_mc   <- mean(full_pop_cvt)
  
  mc_results <- list()
  for (n in this_campaign_grid) {
    cvt_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      loc_idx      <- sample(n_locations, 1)   # pick one location's realized time series
      idx          <- sample(m, n)             # subsample n campaigns from it
      vals         <- sim_mat_time[idx, loc_idx]
      cvt_draws[i] <- sd(vals) / mean(vals)
    }
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n, CVt = cvt_draws,
                                     Ref_CVt_Real = ref_cvt_real, Ref_CVt_MC = ref_cvt_mc)
  }
  temporal_mc_results_list[[cc]] <- bind_rows(mc_results)
}
temporal_mc_results <- bind_rows(temporal_mc_results_list)

temporal_mc_summary <- temporal_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVt = median(CVt), P05_CVt = quantile(CVt, 0.05), P95_CVt = quantile(CVt, 0.95),
            SD_CVt = sd(CVt), Ref_CVt_Real = first(Ref_CVt_Real), Ref_CVt_MC = first(Ref_CVt_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC     = (Median_CVt - Ref_CVt_MC) / Ref_CVt_MC * 100,
         Pct_Bias_vs_Real   = (Median_CVt - Ref_CVt_Real) / Ref_CVt_Real * 100,
         Pct_Gap_Real_vs_MC = (Ref_CVt_Real - Ref_CVt_MC) / Ref_CVt_MC * 100)
write_csv(temporal_mc_summary, file.path(data_out_dir, "CVt_MC_summary_nm.csv"))

bias_annot <- temporal_mc_summary %>% filter(N == actual_n_campaigns) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias (vs model truth)", N, Pct_Bias_vs_MC))

p_cvt <- ggplot(temporal_mc_summary, aes(x = N, y = Median_CVt)) +
  geom_ribbon(aes(ymin = P05_CVt, ymax = P95_CVt), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVt_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_CVt_Real, color = "Realized truth (real sites)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_campaigns, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (real sites)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N campaigns)", y = "CVt (median, 5th-95th pct band)",
       title = "Temporal sample-size sensitivity -- nm",
       subtitle = "Dashed = reference CVt; dotted red = actual sample size") +
  theme_bw() +
  theme(legend.position = "bottom")
p_cvt
ggsave(file.path(plot_dir, "CVt_MC_plot_nm.png"), p_cvt, width = 7, height = 5, dpi = 150)


# compute optimums
bias_tol       <- 10
elbow_fraction <- 0.05

optimal_n_list <- list()
for (cc in constituents) {
  cc_summary <- temporal_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(N)
  
  n_vals    <- cc_summary$N
  width_pct <- (cc_summary$P95_CVt - cc_summary$P05_CVt) / cc_summary$Ref_CVt_MC * 100
  
  marginal_rate   <- -diff(width_pct) / diff(n_vals)
  elbow_idx       <- which(marginal_rate <= elbow_fraction * max(marginal_rate))[1]
  n_opt_precision <- if (is.na(elbow_idx)) NA_real_ else n_vals[elbow_idx + 1]
  
  within_tol   <- abs(cc_summary$Pct_Bias_vs_MC) <= bias_tol
  stays_within <- sapply(seq_along(within_tol), function(i) all(within_tol[i:length(within_tol)]))
  bias_idx     <- which(stays_within)[1]
  n_opt_bias   <- if (is.na(bias_idx)) NA_real_ else n_vals[bias_idx]
  
  optimal_n_list[[cc]] <- tibble(Constituent = cc,
                                 N_optimal_precision = n_opt_precision,
                                 N_optimal_bias      = n_opt_bias)
}
optimal_n_campaigns <- bind_rows(optimal_n_list)
print(optimal_n_campaigns)

#### Testing missing-data tolerance ####
## Pools all real (site x campaign) readings across all real sites and
## randomly drops a % of those INDIVIDUAL readings -- rather than whole
## campaigns -- since with only 3 real campaigns per site, dropping at the
## campaign level would only support 3 coarse outcomes (keep 1, 2, or 3).
## CVt is recomputed per site from whichever campaigns remain there in that
## draw; a site left with fewer than 2 readings is excluded from that draw's
## across-site average.

missing_grid <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8)

missing_mc_results_list <- list()
for (cc in constituents) {
  readings <- clean %>%
    select(Site, Value = all_of(cc)) %>%
    filter(!is.na(Value))
  
  site_factor <- as.integer(factor(readings$Site))
  values      <- readings$Value
  n_readings  <- length(values)
  
  ## Baseline: CVt per site using all real readings (no missingness),
  ## averaged across sites with at least 2 readings.
  ref_cvt_missing <- readings %>%
    group_by(Site) %>%
    filter(n() >= 2) %>%
    summarize(CVt = sd(Value) / mean(Value), .groups = "drop") %>%
    pull(CVt) %>% mean()
  
  mc_results <- list()
  for (pct_missing in missing_grid) {
    n_keep    <- max(1, round(n_readings * (1 - pct_missing)))
    cvt_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      kept_idx <- sample(n_readings, n_keep)
      kept_vals_by_site <- split(values[kept_idx], site_factor[kept_idx])
      site_cvt <- sapply(kept_vals_by_site, function(v) {
        if (length(v) >= 2) sd(v) / mean(v) else NA_real_
      })
      cvt_draws[i] <- if (all(is.na(site_cvt))) NA_real_ else mean(site_cvt, na.rm = TRUE)
    }
    mc_results[[paste(pct_missing)]] <- tibble(Constituent = cc, Pct_Missing = pct_missing * 100,
                                               CVt = cvt_draws, Ref_CVt_Missing = ref_cvt_missing)
  }
  missing_mc_results_list[[cc]] <- bind_rows(mc_results)
}
missing_mc_results <- bind_rows(missing_mc_results_list)

missing_mc_summary <- missing_mc_results %>%
  group_by(Constituent, Pct_Missing) %>%
  summarize(Median_CVt = median(CVt, na.rm = TRUE),
            P05_CVt = quantile(CVt, 0.05, na.rm = TRUE),
            P95_CVt = quantile(CVt, 0.95, na.rm = TRUE),
            SD_CVt = sd(CVt, na.rm = TRUE),
            Pct_Undefined = mean(is.na(CVt)) * 100,
            Ref_CVt_Missing = first(Ref_CVt_Missing), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVt - Ref_CVt_Missing) / Ref_CVt_Missing * 100)
write_csv(missing_mc_summary, file.path(data_out_dir, "CVt_missing_data_summary_nm.csv"))

p_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_CVt)) +
  geom_ribbon(aes(ymin = P05_CVt, ymax = P95_CVt), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVt_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "CVt (median, 5th-95th pct band)",
       title = "Missing-data tolerance (temporal) -- nm",
       subtitle = "Blue dashed = no-missing-data baseline (real sites, 3 real campaigns)") +
  theme_bw()
p_missing
ggsave(file.path(plot_dir, "CVt_missing_data_plot_nm.png"), p_missing, width = 7, height = 5, dpi = 150)


## Compute the maximum allowable % missing data such that bias stays below
## a chosen tolerance, interpolated between the tested missing_grid values

bias_tol_missing  <- 10
width_tol_missing <- 20

find_allowable_pct <- function(pct_missing_vals, metric_vals, tol) {
  pct_missing_vals <- c(0, pct_missing_vals)
  metric_vals      <- c(0, metric_vals)
  if (all(metric_vals <= tol)) {
    max(pct_missing_vals)
  } else if (metric_vals[2] > tol) {
    0
  } else {
    approx(x = metric_vals, y = pct_missing_vals, xout = tol)$y
  }
}

allowable_missing_list <- list()
for (cc in constituents) {
  cc_summary <- missing_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(Pct_Missing) %>%
    mutate(Width_pct = (P95_CVt - P05_CVt) / Ref_CVt_Missing * 100)
  
  pct_allowable_bias  <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$Pct_Bias), bias_tol_missing)
  pct_allowable_width <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_pct, width_tol_missing)
  
  allowable_missing_list[[cc]] <- tibble(Constituent = cc,
                                         Allowable_Pct_Missing_Bias  = pct_allowable_bias,
                                         Allowable_Pct_Missing_Width = pct_allowable_width,
                                         Allowable_Pct_Missing       = min(pct_allowable_bias, pct_allowable_width))
}
allowable_missing <- bind_rows(allowable_missing_list)
print(allowable_missing)


#### ----- br ----- ####
#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_iter        <- 10000               # Monte Carlo iterations
campaign_grid <- c(3:30, 40, 50, 60)  # campaign counts to test

actual_n_campaigns <- 3  # real number of campaigns sampled 

set.seed(42)

clean                   <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
synthetic_extended_full <- read_csv(file.path(data_out_dir, "br_synthetic_with_real.csv"), show_col_types = FALSE)

real_predids <- synthetic_extended_full %>%
  filter(Source == "real") %>%
  pull(predID) %>%
  unique()

synthetic_extended_real_locs <- synthetic_extended_full %>%
  filter(predID %in% real_predids)

#### Varying sample size in time ####

temporal_mc_results_list <- list()
for (cc in constituents) {
  ## Population matrix: rows = synthetic campaigns, columns = locations.
  ## Each column is one location's full realized time series -- the unit we
  ## resample a random one of, per iteration, in place of a fresh MC draw.
  sim_mat_time <- synthetic_extended_full %>%
    filter(Constituent == cc) %>%
    select(CampaignNum, predID, Conc) %>%
    pivot_wider(names_from = predID, values_from = Conc) %>%
    arrange(CampaignNum) %>%
    select(-CampaignNum) %>%
    as.matrix()
  
  m           <- nrow(sim_mat_time)  # number of campaigns in the full synthetic time series
  n_locations <- ncol(sim_mat_time)  # number of locations available to draw from
  this_campaign_grid <- sort(unique(c(campaign_grid[campaign_grid <= m], m)))
  
  ## "Realized truth": CVt at the real-equivalent site locations, computed
  ## over all m campaigns. Tied to the specific locations actually sampled.
  ref_cvt_real <- synthetic_extended_real_locs %>%
    filter(Constituent == cc) %>%
    group_by(predID) %>%
    summarize(CVt = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVt) %>% mean()
  
  ## "Model truth": the full-time-population CVt (all m campaigns), averaged
  ## across every location. Internally consistent with the subsampling loop
  ## below, since it's drawn from the same sim_mat_time.
  full_pop_cvt <- apply(sim_mat_time, 2, function(vals) sd(vals) / mean(vals))
  ref_cvt_mc   <- mean(full_pop_cvt)
  
  mc_results <- list()
  for (n in this_campaign_grid) {
    cvt_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      loc_idx      <- sample(n_locations, 1)   # pick one location's realized time series
      idx          <- sample(m, n)             # subsample n campaigns from it
      vals         <- sim_mat_time[idx, loc_idx]
      cvt_draws[i] <- sd(vals) / mean(vals)
    }
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n, CVt = cvt_draws,
                                     Ref_CVt_Real = ref_cvt_real, Ref_CVt_MC = ref_cvt_mc)
  }
  temporal_mc_results_list[[cc]] <- bind_rows(mc_results)
}
temporal_mc_results <- bind_rows(temporal_mc_results_list)

temporal_mc_summary <- temporal_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVt = median(CVt), P05_CVt = quantile(CVt, 0.05), P95_CVt = quantile(CVt, 0.95),
            SD_CVt = sd(CVt), Ref_CVt_Real = first(Ref_CVt_Real), Ref_CVt_MC = first(Ref_CVt_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC     = (Median_CVt - Ref_CVt_MC) / Ref_CVt_MC * 100,
         Pct_Bias_vs_Real   = (Median_CVt - Ref_CVt_Real) / Ref_CVt_Real * 100,
         Pct_Gap_Real_vs_MC = (Ref_CVt_Real - Ref_CVt_MC) / Ref_CVt_MC * 100)
write_csv(temporal_mc_summary, file.path(data_out_dir, "CVt_MC_summary_br.csv"))

bias_annot <- temporal_mc_summary %>% filter(N == actual_n_campaigns) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias (vs model truth)", N, Pct_Bias_vs_MC))

p_cvt <- ggplot(temporal_mc_summary, aes(x = N, y = Median_CVt)) +
  geom_ribbon(aes(ymin = P05_CVt, ymax = P95_CVt), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVt_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_CVt_Real, color = "Realized truth (real sites)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_campaigns, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (real sites)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N campaigns)", y = "CVt (median, 5th-95th pct band)",
       title = "Temporal sample-size sensitivity -- br",
       subtitle = "Dashed = reference CVt; dotted red = actual sample size") +
  theme_bw() +
  theme(legend.position = "bottom")
p_cvt
ggsave(file.path(plot_dir, "CVt_MC_plot_br.png"), p_cvt, width = 7, height = 5, dpi = 150)


# compute optimums
bias_tol       <- 10
elbow_fraction <- 0.05

optimal_n_list <- list()
for (cc in constituents) {
  cc_summary <- temporal_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(N)
  
  n_vals    <- cc_summary$N
  width_pct <- (cc_summary$P95_CVt - cc_summary$P05_CVt) / cc_summary$Ref_CVt_MC * 100
  
  marginal_rate   <- -diff(width_pct) / diff(n_vals)
  elbow_idx       <- which(marginal_rate <= elbow_fraction * max(marginal_rate))[1]
  n_opt_precision <- if (is.na(elbow_idx)) NA_real_ else n_vals[elbow_idx + 1]
  
  within_tol   <- abs(cc_summary$Pct_Bias_vs_MC) <= bias_tol
  stays_within <- sapply(seq_along(within_tol), function(i) all(within_tol[i:length(within_tol)]))
  bias_idx     <- which(stays_within)[1]
  n_opt_bias   <- if (is.na(bias_idx)) NA_real_ else n_vals[bias_idx]
  
  optimal_n_list[[cc]] <- tibble(Constituent = cc,
                                 N_optimal_precision = n_opt_precision,
                                 N_optimal_bias      = n_opt_bias)
}
optimal_n_campaigns <- bind_rows(optimal_n_list)
print(optimal_n_campaigns)

#### Testing missing-data tolerance ####
## Pools all real (site x campaign) readings across all real sites and
## randomly drops a % of those INDIVIDUAL readings -- rather than whole
## campaigns -- since with only 3 real campaigns per site, dropping at the
## campaign level would only support 3 coarse outcomes (keep 1, 2, or 3).
## CVt is recomputed per site from whichever campaigns remain there in that
## draw; a site left with fewer than 2 readings is excluded from that draw's
## across-site average.

missing_grid <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8)

missing_mc_results_list <- list()
for (cc in constituents) {
  readings <- clean %>%
    select(Site, Value = all_of(cc)) %>%
    filter(!is.na(Value))
  
  site_factor <- as.integer(factor(readings$Site))
  values      <- readings$Value
  n_readings  <- length(values)
  
  ## Baseline: CVt per site using all real readings (no missingness),
  ## averaged across sites with at least 2 readings.
  ref_cvt_missing <- readings %>%
    group_by(Site) %>%
    filter(n() >= 2) %>%
    summarize(CVt = sd(Value) / mean(Value), .groups = "drop") %>%
    pull(CVt) %>% mean()
  
  mc_results <- list()
  for (pct_missing in missing_grid) {
    n_keep    <- max(1, round(n_readings * (1 - pct_missing)))
    cvt_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      kept_idx <- sample(n_readings, n_keep)
      kept_vals_by_site <- split(values[kept_idx], site_factor[kept_idx])
      site_cvt <- sapply(kept_vals_by_site, function(v) {
        if (length(v) >= 2) sd(v) / mean(v) else NA_real_
      })
      cvt_draws[i] <- if (all(is.na(site_cvt))) NA_real_ else mean(site_cvt, na.rm = TRUE)
    }
    mc_results[[paste(pct_missing)]] <- tibble(Constituent = cc, Pct_Missing = pct_missing * 100,
                                               CVt = cvt_draws, Ref_CVt_Missing = ref_cvt_missing)
  }
  missing_mc_results_list[[cc]] <- bind_rows(mc_results)
}
missing_mc_results <- bind_rows(missing_mc_results_list)

missing_mc_summary <- missing_mc_results %>%
  group_by(Constituent, Pct_Missing) %>%
  summarize(Median_CVt = median(CVt, na.rm = TRUE),
            P05_CVt = quantile(CVt, 0.05, na.rm = TRUE),
            P95_CVt = quantile(CVt, 0.95, na.rm = TRUE),
            SD_CVt = sd(CVt, na.rm = TRUE),
            Pct_Undefined = mean(is.na(CVt)) * 100,
            Ref_CVt_Missing = first(Ref_CVt_Missing), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVt - Ref_CVt_Missing) / Ref_CVt_Missing * 100)
write_csv(missing_mc_summary, file.path(data_out_dir, "CVt_missing_data_summary_br.csv"))

p_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_CVt)) +
  geom_ribbon(aes(ymin = P05_CVt, ymax = P95_CVt), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVt_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "CVt (median, 5th-95th pct band)",
       title = "Missing-data tolerance (temporal) -- br",
       subtitle = "Blue dashed = no-missing-data baseline (real sites, 3 real campaigns)") +
  theme_bw()
p_missing
ggsave(file.path(plot_dir, "CVt_missing_data_plot_br.png"), p_missing, width = 7, height = 5, dpi = 150)


## Compute the maximum allowable % missing data such that bias stays below
## a chosen tolerance, interpolated between the tested missing_grid values

bias_tol_missing  <- 10
width_tol_missing <- 20

find_allowable_pct <- function(pct_missing_vals, metric_vals, tol) {
  pct_missing_vals <- c(0, pct_missing_vals)
  metric_vals      <- c(0, metric_vals)
  if (all(metric_vals <= tol)) {
    max(pct_missing_vals)
  } else if (metric_vals[2] > tol) {
    0
  } else {
    approx(x = metric_vals, y = pct_missing_vals, xout = tol)$y
  }
}

allowable_missing_list <- list()
for (cc in constituents) {
  cc_summary <- missing_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(Pct_Missing) %>%
    mutate(Width_pct = (P95_CVt - P05_CVt) / Ref_CVt_Missing * 100)
  
  pct_allowable_bias  <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$Pct_Bias), bias_tol_missing)
  pct_allowable_width <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_pct, width_tol_missing)
  
  allowable_missing_list[[cc]] <- tibble(Constituent = cc,
                                         Allowable_Pct_Missing_Bias  = pct_allowable_bias,
                                         Allowable_Pct_Missing_Width = pct_allowable_width,
                                         Allowable_Pct_Missing       = min(pct_allowable_bias, pct_allowable_width))
}
allowable_missing <- bind_rows(allowable_missing_list)
print(allowable_missing)

