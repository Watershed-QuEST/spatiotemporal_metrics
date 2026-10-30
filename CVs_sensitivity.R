#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-10-28

# Requires: 02_build_synthetic_data.R must be run first (produces data/*_clean.csv, data/*_field_setups.rds, and data/*_synthetic_extended.csv).

# This script performs a sensitivity analysis to answer the question: how much does a CVs estimate vary/undershoot depending on how many sites you sample? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possibile fits to a SSN2 spatial covariance model. 
# Ref_CVs (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R, which simulates a synoptic campaign that sampled every 250 m along the network.

# Outputs:
#   data/CVs_MC_summary.csv:  MC summary
#   plots/CVs_MC_plot.png: CVs vs N, with reference line and bias annotation
#   plots/CVs_sensitivity_metrics_plot.png: percent bias and SD vs N, side by side

#### Packages ####
library(tidyverse)
library(patchwork)

#### ----- nm ----- ####
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
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_with_real.csv"), show_col_types = FALSE)
real_campaign_nums <- synthetic_extended %>%
  filter(Source == "real") %>%
  pull(CampaignNum) %>%
  unique()

synthetic_extended <- synthetic_extended %>%
  filter(CampaignNum %in% real_campaign_nums)

real_predids <- synthetic_extended %>%
  filter(Source == "real") %>%
  pull(predID) %>%
  unique()

#### Look at real data ####

campaign_vals <- split(clean$TDN..mg.N.L., clean$CampaignID)
campaign_vals <- lapply(campaign_vals, function(x) x[!is.na(x)])

for (cc in constituents) {
  campaign_vals <- split(clean[[cc]], clean$CampaignID)
  campaign_vals <- lapply(campaign_vals, function(x) x[!is.na(x)])
  
  cat("\n", cc, "\n")
  for (name in names(campaign_vals)) {
    v <- campaign_vals[[name]]
    cat(" ", name, " mean:", round(mean(v), 3),
        " sd:", round(sd(v), 3),
        " range:", round(diff(range(v)), 3),
        " CV:", round(sd(v) / mean(v), 3), "\n")
  }
}

# plot underlying data distribution
box_data <- clean %>%
  pivot_longer(cols = all_of(constituents), names_to = "Constituent", values_to = "Conc") %>%
  filter(!is.na(Conc))

p_box <- ggplot(box_data, aes(x = CampaignID, y = Conc)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 1) +
  labs(x = "Campaign", y = "Concentration", title = "Boxplot") +
  theme_bw()

p_hist <- ggplot(box_data, aes(x = Conc)) +
  geom_histogram(bins = 10, color = "white") +
  facet_grid(Constituent ~ CampaignID, scales = "free") +
  labs(x = "Concentration", y = "Count", title = "Histograms by campaign") +
  theme_bw()

p_dist <- (p_box | p_hist) +
  plot_annotation(title = "Distribution of real concentration data -- nm")
p_dist
ggsave(file.path(plot_dir, "data_distribution_nm.png"), p_dist, width = 11, height = 5, dpi = 150)

#### Varying sample size in space ####

spatial_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  this_site_grid <- sort(unique(c(site_grid[site_grid <= m], m)))  # always test the full spatial population too
  
  ## "Realized truth": the CVs implied by the 3 real-equivalent campaigns,
  ## computed over all m locations. Tied to the specific conditions on those
  ## 3 real sampling dates.
  ref_cvs_real <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ## "Model truth": the full-population CVs (all m locations), averaged across
  ## MC draws. Internally consistent with the subsampling loop below, since
  ## it's drawn from the exact same sim_mat -- so Median_CVs is guaranteed to
  ## converge to this as N -> m.
  full_pop_cvs <- apply(sim_mat, 2, function(vals) sd(vals) / mean(vals))
  ref_cvs_mc   <- mean(full_pop_cvs)
  
  mc_results <- list()
  for (n in this_site_grid) {
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx          <- sample(m, n)
      vals         <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n, CVs = cvs_draws,
                                     Ref_CVs_Real = ref_cvs_real, Ref_CVs_MC = ref_cvs_mc)
  }
  spatial_mc_results_list[[cc]] <- bind_rows(mc_results)
}
spatial_mc_results <- bind_rows(spatial_mc_results_list)

spatial_mc_summary <- spatial_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs_Real = first(Ref_CVs_Real), Ref_CVs_MC = first(Ref_CVs_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC     = (Median_CVs - Ref_CVs_MC) / Ref_CVs_MC * 100,
         Pct_Bias_vs_Real   = (Median_CVs - Ref_CVs_Real) / Ref_CVs_Real * 100,
         Pct_Gap_Real_vs_MC = (Ref_CVs_Real - Ref_CVs_MC) / Ref_CVs_MC * 100)
write_csv(spatial_mc_summary, file.path(data_out_dir, "CVs_MC_summary_nm.csv"))

bias_annot <- spatial_mc_summary %>% filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias (vs model truth)", N, Pct_Bias_vs_MC))

p_cvs <- ggplot(spatial_mc_summary, aes(x = N, y = Median_CVs)) +
  geom_ribbon(aes(ymin = P05_CVs, ymax = P95_CVs), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVs_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_CVs_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites)", y = "CVs (median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity -- nm",
       subtitle = "Dashed = reference CVs; dotted red = actual sample size") +
  theme_bw()+
  theme(legend.position = "bottom")
p_cvs
ggsave(file.path(plot_dir, "CVs_MC_plot_nm.png"), p_cvs, width = 7, height = 5, dpi = 150)


# compute optimums
bias_tol       <- 10    # percent; acceptable bias vs. model truth
elbow_fraction <- 0.05 # elbow = marginal gain has dropped to this fraction of its steepest value

optimal_n_list <- list()
for (cc in constituents) {
  cc_summary <- spatial_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(N)
  
  n_vals    <- cc_summary$N
  width_pct <- (cc_summary$P95_CVs - cc_summary$P05_CVs) / cc_summary$Ref_CVs_MC * 100
  
  ## Criterion 1: precision elbow
  marginal_rate   <- -diff(width_pct) / diff(n_vals)   # improvement in band width per added site
  elbow_idx       <- which(marginal_rate <= elbow_fraction * max(marginal_rate))[1]
  n_opt_precision <- if (is.na(elbow_idx)) NA_real_ else n_vals[elbow_idx + 1]
  
  ## Criterion 2: bias threshold (vs. model truth; swap to Pct_Bias_vs_Real to
  ## use the realized-truth reference instead)
  within_tol   <- abs(cc_summary$Pct_Bias_vs_MC) <= bias_tol
  stays_within <- sapply(seq_along(within_tol), function(i) all(within_tol[i:length(within_tol)]))
  bias_idx     <- which(stays_within)[1]
  n_opt_bias   <- if (is.na(bias_idx)) NA_real_ else n_vals[bias_idx]
  
  optimal_n_list[[cc]] <- tibble(Constituent = cc,
                                 N_optimal_precision = n_opt_precision,
                                 N_optimal_bias      = n_opt_bias)
}
optimal_n <- bind_rows(optimal_n_list)
print(optimal_n)


#### Testing missing-data tolerance ####
## Starting from the real data as actually collected (3 real campaigns, ~actual_n_sites sites), randomly drop data at increasing missing-data rates and see how much the resulting CVs estimate degrades relative to the no-missing-data baseline.

missing_grid <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8)   # fraction of data missing, tested at random

missing_mc_results_list <- list()
for (cc in constituents) {
  campaign_vals <- split(clean[[cc]], clean$CampaignID)
  campaign_vals <- lapply(campaign_vals, function(x) x[!is.na(x)])
  
  ## Baseline: CVs per real campaign using all the real data (no missingness),
  ## averaged across the 3 real campaigns.
  ref_cvs_missing <- mean(sapply(campaign_vals, function(vals) sd(vals) / mean(vals)))
  
  mc_results <- list()
  for (pct_missing in missing_grid) {
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      campaign_cvs <- sapply(campaign_vals, function(vals) {
        n_keep <- max(2, round(length(vals) * (1 - pct_missing)))
        kept   <- sample(vals, n_keep)
        sd(kept) / mean(kept)
      })
      cvs_draws[i] <- mean(campaign_cvs)
    }
    mc_results[[paste(pct_missing)]] <- tibble(Constituent = cc, Pct_Missing = pct_missing * 100,
                                               CVs = cvs_draws, Ref_CVs_Missing = ref_cvs_missing)
  }
  missing_mc_results_list[[cc]] <- bind_rows(mc_results)
}
missing_mc_results <- bind_rows(missing_mc_results_list)

missing_mc_summary <- missing_mc_results %>%
  group_by(Constituent, Pct_Missing) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs_Missing = first(Ref_CVs_Missing), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVs - Ref_CVs_Missing) / Ref_CVs_Missing * 100)
write_csv(missing_mc_summary, file.path(data_out_dir, "CVs_missing_data_summary_nm.csv"))

p_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_CVs)) +
  geom_ribbon(aes(ymin = P05_CVs, ymax = P95_CVs), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVs_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "CVs (median, 5th-95th pct band)",
       title = "Missing-data tolerance -- nm",
       subtitle = "Blue dashed = no-missing-data baseline (3 real campaigns, ~23 sites)") +
  theme_bw()
p_missing
ggsave(file.path(plot_dir, "CVs_missing_data_plot_nm.png"), p_missing, width = 7, height = 5, dpi = 150)


## Compute the maximum allowable % missing data such that bias stays below a chosen tolerance, interpolated between the tested missing_grid values

bias_tol_missing  <- 10   # percent; acceptable bias vs. no-missing-data baseline
width_tol_missing <- 20   # percent; acceptable width of the P05-P95 band, relative to baseline CVs

## Linear-interpolate the missing-data rate at which a metric (bias or band width, both expressed as % of the baseline CVs) first crosses a tolerance, given values at 0% missing (metric = 0 by definition) plus the tested grid.
find_allowable_pct <- function(pct_missing_vals, metric_vals, tol) {
  pct_missing_vals <- c(0, pct_missing_vals)
  metric_vals      <- c(0, metric_vals)
  if (all(metric_vals <= tol)) {
    max(pct_missing_vals)                 # never crosses tolerance in the tested range
  } else if (metric_vals[2] > tol) {
    0                                      # already exceeds tolerance at the smallest tested rate
  } else {
    approx(x = metric_vals, y = pct_missing_vals, xout = tol)$y
  }
}

allowable_missing_list <- list()
for (cc in constituents) {
  cc_summary <- missing_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(Pct_Missing) %>%
    mutate(Width_pct = (P95_CVs - P05_CVs) / Ref_CVs_Missing * 100)
  
  pct_allowable_bias  <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$Pct_Bias), bias_tol_missing)
  pct_allowable_width <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_pct, width_tol_missing)
  
  allowable_missing_list[[cc]] <- tibble(Constituent = cc,
                                         Allowable_Pct_Missing_Bias  = pct_allowable_bias,
                                         Allowable_Pct_Missing_Width = pct_allowable_width,
                                         Allowable_Pct_Missing       = min(pct_allowable_bias, pct_allowable_width))
}
allowable_missing <- bind_rows(allowable_missing_list)
print(allowable_missing)



#### Varying spatial distribution of sites (mainstem vs. tributary) ####
## Holds N fixed at actual_n_sites and instead varies where those sites sit in the network: tributary-dominated, uniform (matching the real network's own mainstem:tributary mix), or mainstem-dominated. 
#"Mainstem-like" is defined by drainage area (Area_m2, from *_predpts_attributes.csv) as >60% of the total watershed area

predpts_attrs      <- read_csv(file.path(data_out_dir, "nm_predpts_attributes.csv"), show_col_types = FALSE)
area_lookup        <- setNames(predpts_attrs$Area_m2, predpts_attrs$predID)

shift_fraction <- 0.99  # how far to move the real network's mainstem fraction toward each extreme

# ---- Explicit, field-relevant area threshold (replaces population median)
#mainstem_threshold <- 2*1000000  # m^2 (.5 km^2) - tune to whatever scale meaningfully separates mainstem-like from tributary-like reaches here
total_watershed_area <- max(predpts_attrs$Area_m2)
mainstem_frac_of_total <- 0.06  # tune
mainstem_threshold <- mainstem_frac_of_total * total_watershed_area

area_lookup <- setNames(predpts_attrs$Area_m2, predpts_attrs$predID)

# Real sites' own mainstem fraction, using the FULL (un-thinned) classification --
# this should reflect what the real network actually looks like, not the thinned pool
real_is_mainstem <- area_lookup[real_predids] >= mainstem_threshold
p_real <- mean(real_is_mainstem)

# ---- Nearest-neighbor spacing, to help choose a minimum separation distance
nn_dist_mat <- as.matrix(dist(predpts_attrs[, c("X", "Y")]))
diag(nn_dist_mat) <- NA
print(summary(apply(nn_dist_mat, 1, min, na.rm = TRUE)))

min_separation_mainstem_m <- 400
min_separation_trib_m     <- 400   # more spacing between tributary points

mainstem_predids_thin <- decluster_points(mainstem_predids_pool, predpts_attrs, min_separation_mainstem_m, seed = 1)
trib_predids_thin     <- decluster_points(trib_predids_pool,     predpts_attrs, min_separation_trib_m,     seed = 2)

# ---- Decluster: thin each group to mutually separated points
decluster_points <- function(predids, coords_lookup, min_dist, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pts <- coords_lookup %>% filter(predID %in% predids)
  keep <- character(0)
  remaining <- pts
  while (nrow(remaining) > 0) {
    pick <- remaining[sample(nrow(remaining), 1), ]
    keep <- c(keep, pick$predID)
    d <- sqrt((remaining$X - pick$X)^2 + (remaining$Y - pick$Y)^2)
    remaining <- remaining[d >= min_dist, ]
  }
  keep
}

is_mainstem_all <- area_lookup >= mainstem_threshold
mainstem_predids_pool <- names(area_lookup)[is_mainstem_all]
trib_predids_pool     <- names(area_lookup)[!is_mainstem_all]

mainstem_predids_thin <- decluster_points(mainstem_predids_pool, predpts_attrs, min_separation_mainstem_m, seed = 1)
trib_predids_thin     <- decluster_points(trib_predids_pool,     predpts_attrs, min_separation_trib_m, seed = 2)

cat("Mainstem pool:", length(mainstem_predids_pool), "-> thinned to", length(mainstem_predids_thin), "\n")
cat("Tributary pool:", length(trib_predids_pool), "-> thinned to", length(trib_predids_thin), "\n")
cat("Need at least", actual_n_sites, "candidates available in the more restrictive pool.\n")

# ---- Simple map of the mainstem/tributary sampling pools
predpts_attrs <- predpts_attrs %>%
  mutate(
    Group = case_when(
      predID %in% mainstem_predids_thin ~ "Mainstem (sampling pool)",
      predID %in% trib_predids_thin     ~ "Tributary (sampling pool)",
      is_mainstem_all[predID]           ~ "Mainstem (thinned out)",
      TRUE                               ~ "Tributary (thinned out)"
    )
  )
p_network_map <- ggplot(predpts_attrs, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 1.8, alpha = 0.8) +
  coord_equal() +
  scale_color_manual(values = c(
    "Mainstem (sampling pool)"  = "red",
    "Tributary (sampling pool)" = "orange",
    "Mainstem (thinned out)"    = "grey",
    "Tributary (thinned out)"   = "#2e86c1"
  )) +
  labs(title = paste0(toupper("nm"), ": mainstem/tributary sampling pools"),
       x = NULL, y = NULL, color = NULL) +
  theme_minimal()
p_network_map

# per-constituent loop
distribution_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  
  ## Area for each of the m predpts
  area_by_pred <- area_lookup[paste0("P", seq_len(m))]
  is_mainstem  <- area_by_pred >= mainstem_threshold
  mainstem_idx <- match(mainstem_predids_thin, paste0("P", seq_len(m)))
  trib_idx     <- match(trib_predids_thin,     paste0("P", seq_len(m)))
  
  ## The real network's own mix, under this same area threshold.
  p_real <- mean(is_mainstem[real_predids])
  
  scenarios <- tibble(
    Scenario = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"),
    Target_mainstem_frac = c(p_real * (1 - shift_fraction), p_real,
                             p_real + shift_fraction * (1 - p_real))
  )
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ref_cvs_real <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  full_pop_cvs <- apply(sim_mat, 2, function(vals) sd(vals) / mean(vals))
  ref_cvs_mc   <- mean(full_pop_cvs)
  
  mc_results <- list()
  for (s in seq_len(nrow(scenarios))) {
    scen_name  <- scenarios$Scenario[s]
    n_mainstem <- min(max(round(scenarios$Target_mainstem_frac[s] * actual_n_sites), 0), actual_n_sites)
    n_trib     <- actual_n_sites - n_mainstem
    
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx          <- c(sample(mainstem_idx, n_mainstem), sample(trib_idx, n_trib))
      vals         <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    mc_results[[scen_name]] <- tibble(Constituent = cc, Scenario = scen_name, CVs = cvs_draws,
                                      Ref_CVs_Real = ref_cvs_real, Ref_CVs_MC = ref_cvs_mc)
  }
  distribution_mc_results_list[[cc]] <- bind_rows(mc_results)
}
distribution_mc_results <- bind_rows(distribution_mc_results_list)

distribution_mc_summary <- distribution_mc_results %>%
  mutate(Scenario = factor(Scenario, levels = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"))) %>%
  group_by(Constituent, Scenario) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs_Real = first(Ref_CVs_Real), Ref_CVs_MC = first(Ref_CVs_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC   = (Median_CVs - Ref_CVs_MC) / Ref_CVs_MC * 100,
         Pct_Bias_vs_Real = (Median_CVs - Ref_CVs_Real) / Ref_CVs_Real * 100)
write_csv(distribution_mc_summary, file.path(data_out_dir, "CVs_distribution_summary_nm.csv"))

p_dist_cvs <- ggplot(distribution_mc_summary, aes(x = Scenario, y = Median_CVs)) +
  geom_pointrange(aes(ymin = P05_CVs, ymax = P95_CVs)) +
  geom_hline(aes(yintercept = Ref_CVs_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_CVs_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "CVs (median, 5th-95th pct band)",
       title = "Spatial sampling-distribution sensitivity -- nm",
       subtitle = paste0("N fixed at actual_n_sites = ", actual_n_sites, "; dashed = reference CVs")) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
p_dist_cvs
ggsave(file.path(plot_dir, "CVs_distribution_plot_nm.png"), p_dist_cvs, width = 7, height = 5, dpi = 150)
#### ----- br ----- ####
#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_iter    <- 10000                                  # Monte Carlo iterations
site_grid <- c(3:50, 75, 100, 150, 300)              # site counts to test

actual_n_sites <- 19  # real number of br sites sampled -- update if that changes

set.seed(42)

clean            <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
field_setups     <- readRDS(file.path(data_out_dir, "br_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "br_synthetic_with_real.csv"), show_col_types = FALSE)
real_campaign_nums <- synthetic_extended %>%
  filter(Source == "real") %>%
  pull(CampaignNum) %>%
  unique()

synthetic_extended <- synthetic_extended %>%
  filter(CampaignNum %in% real_campaign_nums)

real_predids <- synthetic_extended %>%
  filter(Source == "real") %>%
  pull(predID) %>%
  unique()


#### Look at real data ####

campaign_vals <- split(clean$TDN..mg.N.L., clean$CampaignID)
campaign_vals <- lapply(campaign_vals, function(x) x[!is.na(x)])

for (cc in constituents) {
  campaign_vals <- split(clean[[cc]], clean$CampaignID)
  campaign_vals <- lapply(campaign_vals, function(x) x[!is.na(x)])
  
  cat("\n", cc, "\n")
  for (name in names(campaign_vals)) {
    v <- campaign_vals[[name]]
    cat(" ", name, " mean:", round(mean(v), 3),
        " sd:", round(sd(v), 3),
        " range:", round(diff(range(v)), 3),
        " CV:", round(sd(v) / mean(v), 3), "\n")
  }
}

# plot underlying data distribution
box_data <- clean %>%
  pivot_longer(cols = all_of(constituents), names_to = "Constituent", values_to = "Conc") %>%
  filter(!is.na(Conc))

p_box <- ggplot(box_data, aes(x = CampaignID, y = Conc)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 1) +
  labs(x = "Campaign", y = "Concentration", title = "Boxplot") +
  theme_bw()

p_hist <- ggplot(box_data, aes(x = Conc)) +
  geom_histogram(bins = 10, color = "white") +
  facet_grid(Constituent ~ CampaignID, scales = "free") +
  labs(x = "Concentration", y = "Count", title = "Histograms by campaign") +
  theme_bw()

p_dist <- (p_box | p_hist) +
  plot_annotation(title = "Distribution of real concentration data -- br")
p_dist
ggsave(file.path(plot_dir, "data_distribution_br.png"), p_dist, width = 11, height = 5, dpi = 150)

#### Varying sample size in space ####

spatial_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  this_site_grid <- sort(unique(c(site_grid[site_grid <= m], m)))  # always test the full spatial population too
  
  ## "Realized truth": the CVs implied by the 3 real-equivalent campaigns,
  ## computed over all m locations. Tied to the specific conditions on those
  ## 3 real sampling dates.
  ref_cvs_real <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ## "Model truth": the full-population CVs (all m locations), averaged across
  ## MC draws. Internally consistent with the subsampling loop below, since
  ## it's drawn from the exact same sim_mat -- so Median_CVs is guaranteed to
  ## converge to this as N -> m.
  full_pop_cvs <- apply(sim_mat, 2, function(vals) sd(vals) / mean(vals))
  ref_cvs_mc   <- mean(full_pop_cvs)
  
  mc_results <- list()
  for (n in this_site_grid) {
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx          <- sample(m, n)
      vals         <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n, CVs = cvs_draws,
                                     Ref_CVs_Real = ref_cvs_real, Ref_CVs_MC = ref_cvs_mc)
  }
  spatial_mc_results_list[[cc]] <- bind_rows(mc_results)
}
spatial_mc_results <- bind_rows(spatial_mc_results_list)

spatial_mc_summary <- spatial_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs_Real = first(Ref_CVs_Real), Ref_CVs_MC = first(Ref_CVs_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC     = (Median_CVs - Ref_CVs_MC) / Ref_CVs_MC * 100,
         Pct_Bias_vs_Real   = (Median_CVs - Ref_CVs_Real) / Ref_CVs_Real * 100,
         Pct_Gap_Real_vs_MC = (Ref_CVs_Real - Ref_CVs_MC) / Ref_CVs_MC * 100)
write_csv(spatial_mc_summary, file.path(data_out_dir, "CVs_MC_summary_br.csv"))

bias_annot <- spatial_mc_summary %>% filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.1f%% bias (vs model truth)", N, Pct_Bias_vs_MC))

p_cvs <- ggplot(spatial_mc_summary, aes(x = N, y = Median_CVs)) +
  geom_ribbon(aes(ymin = P05_CVs, ymax = P95_CVs), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVs_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_CVs_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites)", y = "CVs (median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity -- br",
       subtitle = "Dashed = reference CVs; dotted red = actual sample size") +
  theme_bw()+
  theme(legend.position = "bottom")
p_cvs
ggsave(file.path(plot_dir, "CVs_MC_plot_br.png"), p_cvs, width = 7, height = 5, dpi = 150)


# compute optimums
bias_tol       <- 10    # percent; acceptable bias vs. model truth
elbow_fraction <- 0.05 # elbow = marginal gain has dropped to this fraction of its steepest value

optimal_n_list <- list()
for (cc in constituents) {
  cc_summary <- spatial_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(N)
  
  n_vals    <- cc_summary$N
  width_pct <- (cc_summary$P95_CVs - cc_summary$P05_CVs) / cc_summary$Ref_CVs_MC * 100
  
  ## Criterion 1: precision elbow
  marginal_rate   <- -diff(width_pct) / diff(n_vals)   # improvement in band width per added site
  elbow_idx       <- which(marginal_rate <= elbow_fraction * max(marginal_rate))[1]
  n_opt_precision <- if (is.na(elbow_idx)) NA_real_ else n_vals[elbow_idx + 1]
  
  ## Criterion 2: bias threshold (vs. model truth; swap to Pct_Bias_vs_Real to
  ## use the realized-truth reference instead)
  within_tol   <- abs(cc_summary$Pct_Bias_vs_MC) <= bias_tol
  stays_within <- sapply(seq_along(within_tol), function(i) all(within_tol[i:length(within_tol)]))
  bias_idx     <- which(stays_within)[1]
  n_opt_bias   <- if (is.na(bias_idx)) NA_real_ else n_vals[bias_idx]
  
  optimal_n_list[[cc]] <- tibble(Constituent = cc,
                                 N_optimal_precision = n_opt_precision,
                                 N_optimal_bias      = n_opt_bias)
}
optimal_n <- bind_rows(optimal_n_list)
print(optimal_n)

#### Testing missing-data tolerance ####
## Starting from the real data as actually collected (3 real campaigns, ~actual_n_sites sites), randomly drop data at increasing missing-data rates and see how much the resulting CVs estimate degrades relative to the no-missing-data baseline.

missing_grid <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8)   # fraction of data missing, tested at random

missing_mc_results_list <- list()
for (cc in constituents) {
  campaign_vals <- split(clean[[cc]], clean$CampaignID)
  campaign_vals <- lapply(campaign_vals, function(x) x[!is.na(x)])
  
  ## Baseline: CVs per real campaign using all the real data (no missingness),
  ## averaged across the 3 real campaigns.
  ref_cvs_missing <- mean(sapply(campaign_vals, function(vals) sd(vals) / mean(vals)))
  
  mc_results <- list()
  for (pct_missing in missing_grid) {
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      campaign_cvs <- sapply(campaign_vals, function(vals) {
        n_keep <- max(2, round(length(vals) * (1 - pct_missing)))
        kept   <- sample(vals, n_keep)
        sd(kept) / mean(kept)
      })
      cvs_draws[i] <- mean(campaign_cvs)
    }
    mc_results[[paste(pct_missing)]] <- tibble(Constituent = cc, Pct_Missing = pct_missing * 100,
                                               CVs = cvs_draws, Ref_CVs_Missing = ref_cvs_missing)
  }
  missing_mc_results_list[[cc]] <- bind_rows(mc_results)
}
missing_mc_results <- bind_rows(missing_mc_results_list)

missing_mc_summary <- missing_mc_results %>%
  group_by(Constituent, Pct_Missing) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs_Missing = first(Ref_CVs_Missing), .groups = "drop") %>%
  mutate(Pct_Bias = (Median_CVs - Ref_CVs_Missing) / Ref_CVs_Missing * 100)
write_csv(missing_mc_summary, file.path(data_out_dir, "CVs_missing_data_summary_br.csv"))

p_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_CVs)) +
  geom_ribbon(aes(ymin = P05_CVs, ymax = P95_CVs), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_CVs_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "CVs (median, 5th-95th pct band)",
       title = "Missing-data tolerance -- br",
       subtitle = "Blue dashed = no-missing-data baseline (3 real campaigns, ~23 sites)") +
  theme_bw()
p_missing
ggsave(file.path(plot_dir, "CVs_missing_data_plot_br.png"), p_missing, width = 7, height = 5, dpi = 150)


## Compute the maximum allowable % missing data such that bias stays below a chosen tolerance, interpolated between the tested missing_grid values

bias_tol_missing  <- 10   # percent; acceptable bias vs. no-missing-data baseline
width_tol_missing <- 20   # percent; acceptable width of the P05-P95 band, relative to baseline CVs

## Linear-interpolate the missing-data rate at which a metric (bias or band
## width, both expressed as % of the baseline CVs) first crosses a tolerance,
## given values at 0% missing (metric = 0 by definition) plus the tested grid.
find_allowable_pct <- function(pct_missing_vals, metric_vals, tol) {
  pct_missing_vals <- c(0, pct_missing_vals)
  metric_vals      <- c(0, metric_vals)
  if (all(metric_vals <= tol)) {
    max(pct_missing_vals)                 # never crosses tolerance in the tested range
  } else if (metric_vals[2] > tol) {
    0                                      # already exceeds tolerance at the smallest tested rate
  } else {
    approx(x = metric_vals, y = pct_missing_vals, xout = tol)$y
  }
}

allowable_missing_list <- list()
for (cc in constituents) {
  cc_summary <- missing_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(Pct_Missing) %>%
    mutate(Width_pct = (P95_CVs - P05_CVs) / Ref_CVs_Missing * 100)
  
  pct_allowable_bias  <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$Pct_Bias), bias_tol_missing)
  pct_allowable_width <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_pct, width_tol_missing)
  
  allowable_missing_list[[cc]] <- tibble(Constituent = cc,
                                         Allowable_Pct_Missing_Bias  = pct_allowable_bias,
                                         Allowable_Pct_Missing_Width = pct_allowable_width,
                                         Allowable_Pct_Missing       = min(pct_allowable_bias, pct_allowable_width))
}
allowable_missing <- bind_rows(allowable_missing_list)
print(allowable_missing)

#### Varying spatial distribution of sites (mainstem vs. tributary) ####
## Holds N fixed at actual_n_sites and instead varies where those sites sit in the network: tributary-dominated, uniform (matching the real network's own mainstem:tributary mix), or mainstem-dominated. 
#"Mainstem-like" is defined by drainage area (Area_m2, from *_predpts_attributes.csv) as >60% of the total watershed area

predpts_attrs      <- read_csv(file.path(data_out_dir, "br_predpts_attributes.csv"), show_col_types = FALSE)
area_lookup        <- setNames(predpts_attrs$Area_m2, predpts_attrs$predID)

shift_fraction <- 0.99  # how far to move the real network's mainstem fraction toward each extreme

# ---- Explicit, field-relevant area threshold (replaces population median)
#mainstem_threshold <- 2*1000000  # m^2 (.5 km^2) - tune to whatever scale meaningfully separates mainstem-like from tributary-like reaches here
total_watershed_area <- max(predpts_attrs$Area_m2)
mainstem_frac_of_total <- 0.06  # tune
mainstem_threshold <- mainstem_frac_of_total * total_watershed_area

area_lookup <- setNames(predpts_attrs$Area_m2, predpts_attrs$predID)

# Real sites' own mainstem fraction, using the FULL (un-thinned) classification --
# this should reflect what the real network actually looks like, not the thinned pool
real_is_mainstem <- area_lookup[real_predids] >= mainstem_threshold
p_real <- mean(real_is_mainstem)

# ---- Nearest-neighbor spacing, to help choose a minimum separation distance
nn_dist_mat <- as.matrix(dist(predpts_attrs[, c("X", "Y")]))
diag(nn_dist_mat) <- NA
print(summary(apply(nn_dist_mat, 1, min, na.rm = TRUE)))

min_separation_mainstem_m <- 400
min_separation_trib_m     <- 400   # more spacing between tributary points

mainstem_predids_thin <- decluster_points(mainstem_predids_pool, predpts_attrs, min_separation_mainstem_m, seed = 1)
trib_predids_thin     <- decluster_points(trib_predids_pool,     predpts_attrs, min_separation_trib_m,     seed = 2)

# ---- Decluster: thin each group to mutually separated points
decluster_points <- function(predids, coords_lookup, min_dist, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pts <- coords_lookup %>% filter(predID %in% predids)
  keep <- character(0)
  remaining <- pts
  while (nrow(remaining) > 0) {
    pick <- remaining[sample(nrow(remaining), 1), ]
    keep <- c(keep, pick$predID)
    d <- sqrt((remaining$X - pick$X)^2 + (remaining$Y - pick$Y)^2)
    remaining <- remaining[d >= min_dist, ]
  }
  keep
}

is_mainstem_all <- area_lookup >= mainstem_threshold
mainstem_predids_pool <- names(area_lookup)[is_mainstem_all]
trib_predids_pool     <- names(area_lookup)[!is_mainstem_all]

mainstem_predids_thin <- decluster_points(mainstem_predids_pool, predpts_attrs, min_separation_mainstem_m, seed = 1)
trib_predids_thin     <- decluster_points(trib_predids_pool,     predpts_attrs, min_separation_trib_m, seed = 2)

cat("Mainstem pool:", length(mainstem_predids_pool), "-> thinned to", length(mainstem_predids_thin), "\n")
cat("Tributary pool:", length(trib_predids_pool), "-> thinned to", length(trib_predids_thin), "\n")
cat("Need at least", actual_n_sites, "candidates available in the more restrictive pool.\n")

# ---- Simple map of the mainstem/tributary sampling pools
predpts_attrs <- predpts_attrs %>%
  mutate(
    Group = case_when(
      predID %in% mainstem_predids_thin ~ "Mainstem (sampling pool)",
      predID %in% trib_predids_thin     ~ "Tributary (sampling pool)",
      is_mainstem_all[predID]           ~ "Mainstem (thinned out)",
      TRUE                               ~ "Tributary (thinned out)"
    )
  )
p_network_map <- ggplot(predpts_attrs, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 1.8, alpha = 0.8) +
  coord_equal() +
  scale_color_manual(values = c(
    "Mainstem (sampling pool)"  = "red",
    "Tributary (sampling pool)" = "orange",
    "Mainstem (thinned out)"    = "grey",
    "Tributary (thinned out)"   = "#2e86c1"
  )) +
  labs(title = paste0(toupper("br"), ": mainstem/tributary sampling pools"),
       x = NULL, y = NULL, color = NULL) +
  theme_minimal()
p_network_map

# per-constituent loop
distribution_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  
  ## Area for each of the m predpts
  area_by_pred <- area_lookup[paste0("P", seq_len(m))]
  is_mainstem  <- area_by_pred >= mainstem_threshold
  mainstem_idx <- match(mainstem_predids_thin, paste0("P", seq_len(m)))
  trib_idx     <- match(trib_predids_thin,     paste0("P", seq_len(m)))
  
  ## The real network's own mix, under this same area threshold.
  p_real <- mean(is_mainstem[real_predids])
  
  scenarios <- tibble(
    Scenario = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"),
    Target_mainstem_frac = c(p_real * (1 - shift_fraction), p_real,
                             p_real + shift_fraction * (1 - p_real))
  )
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ref_cvs_real <- synthetic_extended %>%
    filter(Constituent == cc) %>%
    group_by(CampaignNum) %>%
    summarize(CVs = sd(Conc) / mean(Conc), .groups = "drop") %>%
    pull(CVs) %>% mean()
  
  full_pop_cvs <- apply(sim_mat, 2, function(vals) sd(vals) / mean(vals))
  ref_cvs_mc   <- mean(full_pop_cvs)
  
  mc_results <- list()
  for (s in seq_len(nrow(scenarios))) {
    scen_name  <- scenarios$Scenario[s]
    n_mainstem <- min(max(round(scenarios$Target_mainstem_frac[s] * actual_n_sites), 0), actual_n_sites)
    n_trib     <- actual_n_sites - n_mainstem
    
    cvs_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx          <- c(sample(mainstem_idx, n_mainstem), sample(trib_idx, n_trib))
      vals         <- sim_mat[idx, i]
      cvs_draws[i] <- sd(vals) / mean(vals)
    }
    mc_results[[scen_name]] <- tibble(Constituent = cc, Scenario = scen_name, CVs = cvs_draws,
                                      Ref_CVs_Real = ref_cvs_real, Ref_CVs_MC = ref_cvs_mc)
  }
  distribution_mc_results_list[[cc]] <- bind_rows(mc_results)
}
distribution_mc_results <- bind_rows(distribution_mc_results_list)

distribution_mc_summary <- distribution_mc_results %>%
  mutate(Scenario = factor(Scenario, levels = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"))) %>%
  group_by(Constituent, Scenario) %>%
  summarize(Median_CVs = median(CVs), P05_CVs = quantile(CVs, 0.05), P95_CVs = quantile(CVs, 0.95),
            SD_CVs = sd(CVs), Ref_CVs_Real = first(Ref_CVs_Real), Ref_CVs_MC = first(Ref_CVs_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC   = (Median_CVs - Ref_CVs_MC) / Ref_CVs_MC * 100,
         Pct_Bias_vs_Real = (Median_CVs - Ref_CVs_Real) / Ref_CVs_Real * 100)
write_csv(distribution_mc_summary, file.path(data_out_dir, "CVs_distribution_summary_br.csv"))

p_dist_cvs <- ggplot(distribution_mc_summary, aes(x = Scenario, y = Median_CVs)) +
  geom_pointrange(aes(ymin = P05_CVs, ymax = P95_CVs)) +
  geom_hline(aes(yintercept = Ref_CVs_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_CVs_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "CVs (median, 5th-95th pct band)",
       title = "Spatial sampling-distribution sensitivity -- br",
       subtitle = paste0("N fixed at actual_n_sites = ", actual_n_sites, "; dashed = reference CVs")) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
p_dist_cvs
ggsave(file.path(plot_dir, "CVs_distribution_plot_br.png"), p_dist_cvs, width = 7, height = 5, dpi = 150)
