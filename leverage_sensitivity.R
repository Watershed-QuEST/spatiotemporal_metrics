#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-11-17
# Last update (Person, Date): 

# Requires: 02_build_synthetic_data.R must be run to produce data/[dataset]_clean.csv for each dataset. If these are available, no need to rerun 02_build_synthetic_data.R

# This script performs sensitivity analyses to answer the questions: 
# 1) how much does a leverage estimate vary depending on how many sites you sample? 
# 2) how much does a leverage estimate vary depending on how much data is missing? 
# 3) how much does a leverage estimate vary depending on how sampling sites are distributed on mainstems vs tributaries? 

# All analyses uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possibile fits to a SSN2 spatial covariance model. 

# Outputs:
#   data/leverage_MC_summary.csv:  MC summary
#   plots/leverage_MC_plot.png: CVs vs N, with reference line and bias annotation
#   plots/leverage_sensitivity_metrics_plot.png: percent bias and SD vs N, side by side

#### Packages ####
library(tidyverse)
library(patchwork)
library(watershedmetrics)
# import variance-collapse breakpoint detector (a script helper using the changepoint package)
source(system.file("scripts", "variance_breakpoints.R", package = "watershedmetrics"))
theme_set(theme_classic(base_size = 11))

#### ----- nm ----- ####
#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_iter    <- 10000
site_grid <- c(3:50, 75, 100, 150, 300)

actual_n_sites  <- 23
outlet_site_id  <- "USF12"

set.seed(42)

clean              <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups       <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_with_real.csv"), show_col_types = FALSE)
predpts_attrs      <- read_csv(file.path(data_out_dir, "nm_predpts_attributes.csv"), show_col_types = FALSE)
area_lookup        <- setNames(predpts_attrs$Area_m2, predpts_attrs$predID)
outlet_predID      <- names(area_lookup)[which.max(area_lookup)]

real_campaign_nums <- synthetic_extended %>%
  filter(Source == "real") %>% pull(CampaignNum) %>% unique()

synthetic_extended <- synthetic_extended %>%
  filter(CampaignNum %in% real_campaign_nums)

real_predids <- synthetic_extended %>%
  filter(Source == "real") %>% pull(predID) %>% unique()

## Real-data specific discharge (L/sec per m^2). mean_real_specific_discharge
## is computed for documentation/reference but isn't multiplied into the
## synthetic-network math below -- with a single representative value applied
## everywhere, q_s/q_o = 1 by construction, so the term cancels out.
clean_prepped <- clean %>%
  mutate(Watershed = "nm", SpecificDischarge = Q / Area.m2)

mean_real_specific_discharge <- mean(clean_prepped$SpecificDischarge, na.rm = TRUE)

## Real-data leverage (true Lpct, measured discharge), computed once for all
## constituents via the wide-format interface.
real_leverage <- subcatchment_leverage(
  data = clean_prepped,
  concentration_col = setNames(constituents, constituents),
  site_col = "Site", event_col = "CampaignID",
  watershed_col = "Watershed",
  contributing_area_col = "Area.m2", outlet_site = outlet_site_id,
  variant = "Lpct", specific_discharge_col = "SpecificDischarge"
)

## "Realized truth": mean (and mean |leverage|) across non-outlet real sites,
## averaged within each real campaign first, then across campaigns.
ref_leverage_real_by_cc <- real_leverage %>%
  filter(Site != outlet_site_id) %>%
  group_by(Constituent, CampaignID) %>%
  summarize(Leverage = mean(leverage, na.rm = TRUE),
            AbsLeverage = mean(abs(leverage), na.rm = TRUE), .groups = "drop") %>%
  group_by(Constituent) %>%
  summarize(Ref_Leverage_Real = mean(Leverage),
            Ref_AbsLeverage_Real = mean(AbsLeverage), .groups = "drop")


#### Varying sample size in space ####

spatial_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  this_site_grid <- sort(unique(c(site_grid[site_grid <= m], m)))
  
  area_by_pred <- area_lookup[paste0("P", seq_len(m))]
  outlet_idx   <- which.max(area_by_pred)
  
  ref_leverage_real     <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_Leverage_Real)
  ref_abs_leverage_real <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_AbsLeverage_Real)
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ## "Model truth": full-population leverage (every non-outlet location vs.
  ## the fixed outlet), averaged across MC draws -- signed and |leverage|.
  full_pop_leverage <- apply(sim_mat, 2, function(vals) {
    lev_pct <- 100 * (vals - vals[outlet_idx]) * (area_by_pred / area_by_pred[outlet_idx]) / vals[outlet_idx]
    c(Leverage = mean(lev_pct[-outlet_idx]), AbsLeverage = mean(abs(lev_pct[-outlet_idx])))
  })
  ref_leverage_mc     <- mean(full_pop_leverage["Leverage", ])
  ref_abs_leverage_mc <- mean(full_pop_leverage["AbsLeverage", ])
  
  non_outlet_pool <- setdiff(seq_len(m), outlet_idx)
  
  mc_results <- list()
  for (n in this_site_grid) {
    lev_draws     <- numeric(n_iter)
    abs_lev_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx     <- c(outlet_idx, sample(non_outlet_pool, n - 1))
      vals    <- sim_mat[idx, i]
      areas   <- area_by_pred[idx]
      lev_pct <- 100 * (vals - vals[1]) * (areas / areas[1]) / vals[1]
      lev_draws[i]     <- mean(lev_pct[-1])
      abs_lev_draws[i] <- mean(abs(lev_pct[-1]))
    }
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n,
                                     Leverage = lev_draws, AbsLeverage = abs_lev_draws,
                                     Ref_Leverage_Real = ref_leverage_real, Ref_Leverage_MC = ref_leverage_mc,
                                     Ref_AbsLeverage_Real = ref_abs_leverage_real, Ref_AbsLeverage_MC = ref_abs_leverage_mc)
  }
  spatial_mc_results_list[[cc]] <- bind_rows(mc_results)
}
spatial_mc_results <- bind_rows(spatial_mc_results_list)

spatial_mc_summary <- spatial_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_Leverage = median(Leverage), P05_Leverage = quantile(Leverage, 0.05),
            P95_Leverage = quantile(Leverage, 0.95), SD_Leverage = sd(Leverage),
            Median_AbsLeverage = median(AbsLeverage), P05_AbsLeverage = quantile(AbsLeverage, 0.05),
            P95_AbsLeverage = quantile(AbsLeverage, 0.95), SD_AbsLeverage = sd(AbsLeverage),
            Ref_Leverage_Real = first(Ref_Leverage_Real), Ref_Leverage_MC = first(Ref_Leverage_MC),
            Ref_AbsLeverage_Real = first(Ref_AbsLeverage_Real), Ref_AbsLeverage_MC = first(Ref_AbsLeverage_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC     = (Median_Leverage - Ref_Leverage_MC) / abs(Ref_Leverage_MC) * 100,
         Abs_Bias_vs_MC     = Median_Leverage - Ref_Leverage_MC,
         Pct_Bias_vs_Real   = (Median_Leverage - Ref_Leverage_Real) / abs(Ref_Leverage_Real) * 100,
         Pct_Gap_Real_vs_MC = (Ref_Leverage_Real - Ref_Leverage_MC) / abs(Ref_Leverage_MC) * 100,
         AbsPct_Bias_vs_MC  = (Median_AbsLeverage - Ref_AbsLeverage_MC) / abs(Ref_AbsLeverage_MC) * 100,
         AbsAbs_Bias_vs_MC  = Median_AbsLeverage - Ref_AbsLeverage_MC)
write_csv(spatial_mc_summary, file.path(data_out_dir, "leverage_MC_summary_nm.csv"))

bias_annot <- spatial_mc_summary %>% filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.2f pct-pts (%+.1f%%) bias", N, Abs_Bias_vs_MC, Pct_Bias_vs_MC))

p_leverage <- ggplot(spatial_mc_summary, aes(x = N, y = Median_Leverage)) +
  geom_ribbon(aes(ymin = P05_Leverage, ymax = P95_Leverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_Leverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_Leverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites, outlet included)", y = "Leverage (%, median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity (leverage) -- nm",
       subtitle = "Dashed = reference leverage; dotted red = actual sample size") +
  theme_bw() +
  theme(legend.position = "bottom")
p_leverage
ggsave(file.path(plot_dir, "leverage_MC_plot_nm.png"), p_leverage, width = 7, height = 5, dpi = 150)

bias_annot_abs <- spatial_mc_summary %>% filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.2f pct-pts (%+.1f%%) bias", N, AbsAbs_Bias_vs_MC, AbsPct_Bias_vs_MC))

p_abs_leverage <- ggplot(spatial_mc_summary, aes(x = N, y = Median_AbsLeverage)) +
  geom_ribbon(aes(ymin = P05_AbsLeverage, ymax = P95_AbsLeverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_AbsLeverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_AbsLeverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot_abs, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites, outlet included)", y = "|Leverage| (%, median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity (|leverage|) -- nm",
       subtitle = "Dashed = reference |leverage|; dotted red = actual sample size") +
  theme_bw() +
  theme(legend.position = "bottom")
p_abs_leverage
ggsave(file.path(plot_dir, "abs_leverage_MC_plot_nm.png"), p_abs_leverage, width = 7, height = 5, dpi = 150)


# compute optimums
bias_tol_abs   <- 5    # allowable |bias| in leverage percentage-points (vs. model truth) -- tune to taste
elbow_fraction <- 0.05

optimal_n_list <- list()
for (cc in constituents) {
  cc_summary <- spatial_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(N)
  
  n_vals        <- cc_summary$N
  width_pct     <- (cc_summary$P95_Leverage - cc_summary$P05_Leverage) / abs(cc_summary$Ref_Leverage_MC) * 100
  width_pct_abs <- (cc_summary$P95_AbsLeverage - cc_summary$P05_AbsLeverage) / abs(cc_summary$Ref_AbsLeverage_MC) * 100
  
  marginal_rate   <- -diff(width_pct) / diff(n_vals)
  elbow_idx       <- which(marginal_rate <= elbow_fraction * max(marginal_rate))[1]
  n_opt_precision <- if (is.na(elbow_idx)) NA_real_ else n_vals[elbow_idx + 1]
  
  marginal_rate_abs   <- -diff(width_pct_abs) / diff(n_vals)
  elbow_idx_abs       <- which(marginal_rate_abs <= elbow_fraction * max(marginal_rate_abs))[1]
  n_opt_precision_abs <- if (is.na(elbow_idx_abs)) NA_real_ else n_vals[elbow_idx_abs + 1]
  
  within_tol   <- abs(cc_summary$Abs_Bias_vs_MC) <= bias_tol_abs
  stays_within <- sapply(seq_along(within_tol), function(i) all(within_tol[i:length(within_tol)]))
  bias_idx     <- which(stays_within)[1]
  n_opt_bias   <- if (is.na(bias_idx)) NA_real_ else n_vals[bias_idx]
  
  within_tol_abs   <- abs(cc_summary$AbsAbs_Bias_vs_MC) <= bias_tol_abs
  stays_within_abs <- sapply(seq_along(within_tol_abs), function(i) all(within_tol_abs[i:length(within_tol_abs)]))
  bias_idx_abs     <- which(stays_within_abs)[1]
  n_opt_bias_abs   <- if (is.na(bias_idx_abs)) NA_real_ else n_vals[bias_idx_abs]
  
  optimal_n_list[[cc]] <- tibble(Constituent = cc,
                                 N_optimal_precision     = n_opt_precision,
                                 N_optimal_precision_abs = n_opt_precision_abs,
                                 N_optimal_bias          = n_opt_bias,
                                 N_optimal_bias_abs      = n_opt_bias_abs)
}
optimal_n <- bind_rows(optimal_n_list)
print(optimal_n)


#### Testing missing-data tolerance ####

missing_grid <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8)

missing_mc_results_list <- list()
for (cc in constituents) {
  cc_data <- clean_prepped %>%
    filter(!is.na(.data[[cc]])) %>%
    select(Site, CampaignID, Conc = all_of(cc), Area.m2, SpecificDischarge)
  
  campaign_outlet <- split(cc_data %>% filter(Site == outlet_site_id),
                           (cc_data %>% filter(Site == outlet_site_id))$CampaignID)
  campaign_sites  <- split(cc_data %>% filter(Site != outlet_site_id),
                           (cc_data %>% filter(Site != outlet_site_id))$CampaignID)
  
  ref_campaign_lev <- sapply(names(campaign_sites), function(camp) {
    sites  <- campaign_sites[[camp]]
    outlet <- campaign_outlet[[camp]]
    lev <- 100 * (sites$Conc - outlet$Conc) * (sites$Area.m2 / outlet$Area.m2) *
      (sites$SpecificDischarge / outlet$SpecificDischarge) / outlet$Conc
    c(Leverage = mean(lev), AbsLeverage = mean(abs(lev)))
  })
  ref_leverage_missing     <- mean(ref_campaign_lev["Leverage", ])
  ref_abs_leverage_missing <- mean(ref_campaign_lev["AbsLeverage", ])
  
  mc_results <- list()
  for (pct_missing in missing_grid) {
    lev_draws     <- numeric(n_iter)
    abs_lev_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      campaign_lev <- sapply(names(campaign_sites), function(camp) {
        sites  <- campaign_sites[[camp]]
        outlet <- campaign_outlet[[camp]]
        n_keep <- max(1, round(nrow(sites) * (1 - pct_missing)))
        kept   <- sites[sample(nrow(sites), n_keep), ]
        lev <- 100 * (kept$Conc - outlet$Conc) * (kept$Area.m2 / outlet$Area.m2) *
          (kept$SpecificDischarge / outlet$SpecificDischarge) / outlet$Conc
        c(Leverage = mean(lev), AbsLeverage = mean(abs(lev)))
      })
      lev_draws[i]     <- mean(campaign_lev["Leverage", ])
      abs_lev_draws[i] <- mean(campaign_lev["AbsLeverage", ])
    }
    mc_results[[paste(pct_missing)]] <- tibble(Constituent = cc, Pct_Missing = pct_missing * 100,
                                               Leverage = lev_draws, AbsLeverage = abs_lev_draws,
                                               Ref_Leverage_Missing = ref_leverage_missing,
                                               Ref_AbsLeverage_Missing = ref_abs_leverage_missing)
  }
  missing_mc_results_list[[cc]] <- bind_rows(mc_results)
}
missing_mc_results <- bind_rows(missing_mc_results_list)

missing_mc_summary <- missing_mc_results %>%
  group_by(Constituent, Pct_Missing) %>%
  summarize(Median_Leverage = median(Leverage), P05_Leverage = quantile(Leverage, 0.05),
            P95_Leverage = quantile(Leverage, 0.95), SD_Leverage = sd(Leverage),
            Median_AbsLeverage = median(AbsLeverage), P05_AbsLeverage = quantile(AbsLeverage, 0.05),
            P95_AbsLeverage = quantile(AbsLeverage, 0.95), SD_AbsLeverage = sd(AbsLeverage),
            Ref_Leverage_Missing = first(Ref_Leverage_Missing),
            Ref_AbsLeverage_Missing = first(Ref_AbsLeverage_Missing), .groups = "drop") %>%
  mutate(Pct_Bias    = (Median_Leverage - Ref_Leverage_Missing) / abs(Ref_Leverage_Missing) * 100,
         Abs_Bias    = Median_Leverage - Ref_Leverage_Missing,
         AbsPct_Bias = (Median_AbsLeverage - Ref_AbsLeverage_Missing) / abs(Ref_AbsLeverage_Missing) * 100,
         AbsAbs_Bias = Median_AbsLeverage - Ref_AbsLeverage_Missing)
write_csv(missing_mc_summary, file.path(data_out_dir, "leverage_missing_data_summary_nm.csv"))

p_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_Leverage)) +
  geom_ribbon(aes(ymin = P05_Leverage, ymax = P95_Leverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_Leverage_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "Leverage (%, median, 5th-95th pct band)",
       title = "Missing-data tolerance (leverage) -- nm",
       subtitle = "Blue dashed = no-missing-data baseline (3 real campaigns); outlet never dropped") +
  theme_bw()
p_missing
ggsave(file.path(plot_dir, "leverage_missing_data_plot_nm.png"), p_missing, width = 7, height = 5, dpi = 150)

p_abs_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_AbsLeverage)) +
  geom_ribbon(aes(ymin = P05_AbsLeverage, ymax = P95_AbsLeverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_AbsLeverage_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "|Leverage| (%, median, 5th-95th pct band)",
       title = "Missing-data tolerance (|leverage|) -- nm",
       subtitle = "Blue dashed = no-missing-data baseline (3 real campaigns); outlet never dropped") +
  theme_bw()
p_abs_missing
ggsave(file.path(plot_dir, "abs_leverage_missing_data_plot_nm.png"), p_abs_missing, width = 7, height = 5, dpi = 150)

### calculate optimums 

bias_tol_missing_abs  <- 1      # leverage percentage-point bias tolerance -- tune to taste
width_tol_missing_abs <- 5      # leverage percentage-point width tolerance -- tune to taste

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
    mutate(Width_abs     = P95_Leverage - P05_Leverage,
           Width_abs_Lev = P95_AbsLeverage - P05_AbsLeverage)
  
  pct_allowable_bias         <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$Abs_Bias), bias_tol_missing_abs)
  pct_allowable_width        <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_abs, width_tol_missing_abs)
  pct_allowable_bias_absLev  <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$AbsAbs_Bias), bias_tol_missing_abs)
  pct_allowable_width_absLev <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_abs_Lev, width_tol_missing_abs)
  
  allowable_missing_list[[cc]] <- tibble(Constituent = cc,
                                         Allowable_Pct_Missing_Bias         = pct_allowable_bias,
                                         Allowable_Pct_Missing_Width        = pct_allowable_width,
                                         #Allowable_Pct_Missing              = min(pct_allowable_bias, pct_allowable_width),
                                         Allowable_Pct_Missing_Bias_AbsLev  = pct_allowable_bias_absLev,
                                         Allowable_Pct_Missing_Width_AbsLev = pct_allowable_width_absLev)
                                         #Allowable_Pct_Missing_AbsLev       = min(pct_allowable_bias_absLev, pct_allowable_width_absLev))
}
allowable_missing <- bind_rows(allowable_missing_list)
print(allowable_missing)

#
#### Varying spatial distribution of sites (mainstem vs. tributary) ####

shift_fraction <- 0.99

total_watershed_area  <- max(predpts_attrs$Area_m2)
mainstem_frac_of_total <- 0.06
mainstem_threshold     <- mainstem_frac_of_total * total_watershed_area

real_is_mainstem <- area_lookup[real_predids] >= mainstem_threshold
p_real <- mean(real_is_mainstem)

nn_dist_mat <- as.matrix(dist(predpts_attrs[, c("X", "Y")]))
diag(nn_dist_mat) <- NA
print(summary(apply(nn_dist_mat, 1, min, na.rm = TRUE)))

min_separation_mainstem_m <- 350
min_separation_trib_m     <- 400

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
mainstem_predids_pool <- setdiff(names(area_lookup)[is_mainstem_all], outlet_predID)
trib_predids_pool     <- setdiff(names(area_lookup)[!is_mainstem_all], outlet_predID)

mainstem_predids_thin <- decluster_points(mainstem_predids_pool, predpts_attrs, min_separation_mainstem_m, seed = 1)
trib_predids_thin     <- decluster_points(trib_predids_pool,     predpts_attrs, min_separation_trib_m, seed = 2)

cat("Mainstem pool:", length(mainstem_predids_pool), "-> thinned to", length(mainstem_predids_thin), "\n")
cat("Tributary pool:", length(trib_predids_pool), "-> thinned to", length(trib_predids_thin), "\n")
cat("Need at least", actual_n_sites - 1, "non-outlet candidates available in the more restrictive pool.\n")

predpts_attrs <- predpts_attrs %>%
  mutate(
    Group = case_when(
      predID == outlet_predID           ~ "Outlet",
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
    "Outlet"                     = "black",
    "Mainstem (sampling pool)"  = "red",
    "Tributary (sampling pool)" = "orange",
    "Mainstem (thinned out)"    = "grey",
    "Tributary (thinned out)"   = "#2e86c1"
  )) +
  labs(title = paste0(toupper("nm"), ": mainstem/tributary sampling pools"),
       x = NULL, y = NULL, color = NULL) +
  theme_minimal()
p_network_map

distribution_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  
  area_by_pred <- area_lookup[paste0("P", seq_len(m))]
  is_mainstem  <- area_by_pred >= mainstem_threshold
  outlet_idx   <- which(paste0("P", seq_len(m)) == outlet_predID)
  mainstem_idx <- match(mainstem_predids_thin, paste0("P", seq_len(m)))
  trib_idx     <- match(trib_predids_thin,     paste0("P", seq_len(m)))
  
  p_real <- mean(is_mainstem[real_predids])
  
  scenarios <- tibble(
    Scenario = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"),
    Target_mainstem_frac = c(p_real * (1 - shift_fraction), p_real,
                             p_real + shift_fraction * (1 - p_real))
  )
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ref_leverage_real     <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_Leverage_Real)
  ref_abs_leverage_real <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_AbsLeverage_Real)
  
  full_pop_leverage <- apply(sim_mat, 2, function(vals) {
    lev_pct <- 100 * (vals - vals[outlet_idx]) * (area_by_pred / area_by_pred[outlet_idx]) / vals[outlet_idx]
    c(Leverage = mean(lev_pct[-outlet_idx]), AbsLeverage = mean(abs(lev_pct[-outlet_idx])))
  })
  ref_leverage_mc     <- mean(full_pop_leverage["Leverage", ])
  ref_abs_leverage_mc <- mean(full_pop_leverage["AbsLeverage", ])
  
  mc_results <- list()
  for (s in seq_len(nrow(scenarios))) {
    scen_name  <- scenarios$Scenario[s]
    n_mainstem <- min(max(round(scenarios$Target_mainstem_frac[s] * (actual_n_sites - 1)), 0), actual_n_sites - 1)
    n_trib     <- (actual_n_sites - 1) - n_mainstem
    
    lev_draws     <- numeric(n_iter)
    abs_lev_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx     <- c(outlet_idx, sample(mainstem_idx, n_mainstem), sample(trib_idx, n_trib))
      vals    <- sim_mat[idx, i]
      areas   <- area_by_pred[idx]
      lev_pct <- 100 * (vals - vals[1]) * (areas / areas[1]) / vals[1]
      lev_draws[i]     <- mean(lev_pct[-1])
      abs_lev_draws[i] <- mean(abs(lev_pct[-1]))
    }
    mc_results[[scen_name]] <- tibble(Constituent = cc, Scenario = scen_name,
                                      Leverage = lev_draws, AbsLeverage = abs_lev_draws,
                                      Ref_Leverage_Real = ref_leverage_real, Ref_Leverage_MC = ref_leverage_mc,
                                      Ref_AbsLeverage_Real = ref_abs_leverage_real, Ref_AbsLeverage_MC = ref_abs_leverage_mc)
  }
  distribution_mc_results_list[[cc]] <- bind_rows(mc_results)
}
distribution_mc_results <- bind_rows(distribution_mc_results_list)

distribution_mc_summary <- distribution_mc_results %>%
  mutate(Scenario = factor(Scenario, levels = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"))) %>%
  group_by(Constituent, Scenario) %>%
  summarize(Median_Leverage = median(Leverage), P05_Leverage = quantile(Leverage, 0.05),
            P95_Leverage = quantile(Leverage, 0.95), SD_Leverage = sd(Leverage),
            Median_AbsLeverage = median(AbsLeverage), P05_AbsLeverage = quantile(AbsLeverage, 0.05),
            P95_AbsLeverage = quantile(AbsLeverage, 0.95), SD_AbsLeverage = sd(AbsLeverage),
            Ref_Leverage_Real = first(Ref_Leverage_Real), Ref_Leverage_MC = first(Ref_Leverage_MC),
            Ref_AbsLeverage_Real = first(Ref_AbsLeverage_Real), Ref_AbsLeverage_MC = first(Ref_AbsLeverage_MC),
            .groups = "drop") %>%
  mutate(Abs_Bias_vs_MC        = Median_Leverage - Ref_Leverage_MC,
         Abs_Bias_vs_Real      = Median_Leverage - Ref_Leverage_Real,
         AbsAbs_Bias_vs_MC     = Median_AbsLeverage - Ref_AbsLeverage_MC,
         AbsAbs_Bias_vs_Real   = Median_AbsLeverage - Ref_AbsLeverage_Real,
         Width_abs             = P95_Leverage - P05_Leverage,
         Width_abs_AbsLeverage = P95_AbsLeverage - P05_AbsLeverage)
write_csv(distribution_mc_summary, file.path(data_out_dir, "leverage_distribution_summary_nm.csv"))

p_dist_leverage <- ggplot(distribution_mc_summary, aes(x = Scenario, y = Median_Leverage)) +
  geom_pointrange(aes(ymin = P05_Leverage, ymax = P95_Leverage)) +
  geom_hline(aes(yintercept = Ref_Leverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_Leverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Leverage (%, median, 5th-95th pct band)",
       title = "Spatial sampling-distribution sensitivity (leverage) -- nm",
       subtitle = paste0("N fixed at actual_n_sites = ", actual_n_sites, " (outlet included); dashed = reference leverage")) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
p_dist_leverage
ggsave(file.path(plot_dir, "leverage_distribution_plot_nm.png"), p_dist_leverage, width = 7, height = 5, dpi = 150)

p_abs_dist_leverage <- ggplot(distribution_mc_summary, aes(x = Scenario, y = Median_AbsLeverage)) +
  geom_pointrange(aes(ymin = P05_AbsLeverage, ymax = P95_AbsLeverage)) +
  geom_hline(aes(yintercept = Ref_AbsLeverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_AbsLeverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "|Leverage| (%, median, 5th-95th pct band)",
       title = "Spatial sampling-distribution sensitivity (|leverage|) -- nm",
       subtitle = paste0("N fixed at actual_n_sites = ", actual_n_sites, " (outlet included); dashed = reference |leverage|")) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
p_abs_dist_leverage
ggsave(file.path(plot_dir, "abs_leverage_distribution_plot_nm.png"), p_abs_dist_leverage, width = 7, height = 5, dpi = 150)

#
#### ----- br ----- ####
#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

n_iter    <- 10000
site_grid <- c(3:50, 75, 100, 150, 300)

actual_n_sites  <- 19
outlet_site_id  <- "BRM01"

set.seed(42)

clean              <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
field_setups       <- readRDS(file.path(data_out_dir, "br_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "br_synthetic_with_real.csv"), show_col_types = FALSE)
predpts_attrs      <- read_csv(file.path(data_out_dir, "br_predpts_attributes.csv"), show_col_types = FALSE)
area_lookup        <- setNames(predpts_attrs$Area_m2, predpts_attrs$predID)
outlet_predID      <- names(area_lookup)[which.max(area_lookup)]

real_campaign_nums <- synthetic_extended %>%
  filter(Source == "real") %>% pull(CampaignNum) %>% unique()

synthetic_extended <- synthetic_extended %>%
  filter(CampaignNum %in% real_campaign_nums)

real_predids <- synthetic_extended %>%
  filter(Source == "real") %>% pull(predID) %>% unique()

## Real-data specific discharge (L/sec per m^2). mean_real_specific_discharge
## is computed for documentation/reference but isn't multiplied into the
## synthetic-network math below -- with a single representative value applied
## everywhere, q_s/q_o = 1 by construction, so the term cancels out.
clean_prepped <- clean %>%
  mutate(Watershed = "br", SpecificDischarge = Q / Area.m2)

mean_real_specific_discharge <- mean(clean_prepped$SpecificDischarge, na.rm = TRUE)

## Real-data leverage (true Lpct, measured discharge), computed once for all
## constituents via the wide-format interface.
real_leverage <- subcatchment_leverage(
  data = clean_prepped,
  concentration_col = setNames(constituents, constituents),
  site_col = "Site", event_col = "CampaignID",
  watershed_col = "Watershed",
  contributing_area_col = "Area.m2", outlet_site = outlet_site_id,
  variant = "Lpct", specific_discharge_col = "SpecificDischarge"
)

## "Realized truth": mean (and mean |leverage|) across non-outlet real sites,
## averaged within each real campaign first, then across campaigns.
ref_leverage_real_by_cc <- real_leverage %>%
  filter(Site != outlet_site_id) %>%
  group_by(Constituent, CampaignID) %>%
  summarize(Leverage = mean(leverage, na.rm = TRUE),
            AbsLeverage = mean(abs(leverage), na.rm = TRUE), .groups = "drop") %>%
  group_by(Constituent) %>%
  summarize(Ref_Leverage_Real = mean(Leverage),
            Ref_AbsLeverage_Real = mean(AbsLeverage), .groups = "drop")


#### Varying sample size in space ####

spatial_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  this_site_grid <- sort(unique(c(site_grid[site_grid <= m], m)))
  
  area_by_pred <- area_lookup[paste0("P", seq_len(m))]
  outlet_idx   <- which.max(area_by_pred)
  
  ref_leverage_real     <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_Leverage_Real)
  ref_abs_leverage_real <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_AbsLeverage_Real)
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ## "Model truth": full-population leverage (every non-outlet location vs.
  ## the fixed outlet), averaged across MC draws -- signed and |leverage|.
  full_pop_leverage <- apply(sim_mat, 2, function(vals) {
    lev_pct <- 100 * (vals - vals[outlet_idx]) * (area_by_pred / area_by_pred[outlet_idx]) / vals[outlet_idx]
    c(Leverage = mean(lev_pct[-outlet_idx]), AbsLeverage = mean(abs(lev_pct[-outlet_idx])))
  })
  ref_leverage_mc     <- mean(full_pop_leverage["Leverage", ])
  ref_abs_leverage_mc <- mean(full_pop_leverage["AbsLeverage", ])
  
  non_outlet_pool <- setdiff(seq_len(m), outlet_idx)
  
  mc_results <- list()
  for (n in this_site_grid) {
    lev_draws     <- numeric(n_iter)
    abs_lev_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx     <- c(outlet_idx, sample(non_outlet_pool, n - 1))
      vals    <- sim_mat[idx, i]
      areas   <- area_by_pred[idx]
      lev_pct <- 100 * (vals - vals[1]) * (areas / areas[1]) / vals[1]
      lev_draws[i]     <- mean(lev_pct[-1])
      abs_lev_draws[i] <- mean(abs(lev_pct[-1]))
    }
    mc_results[[paste(n)]] <- tibble(Constituent = cc, N = n,
                                     Leverage = lev_draws, AbsLeverage = abs_lev_draws,
                                     Ref_Leverage_Real = ref_leverage_real, Ref_Leverage_MC = ref_leverage_mc,
                                     Ref_AbsLeverage_Real = ref_abs_leverage_real, Ref_AbsLeverage_MC = ref_abs_leverage_mc)
  }
  spatial_mc_results_list[[cc]] <- bind_rows(mc_results)
}
spatial_mc_results <- bind_rows(spatial_mc_results_list)

spatial_mc_summary <- spatial_mc_results %>%
  group_by(Constituent, N) %>%
  summarize(Median_Leverage = median(Leverage), P05_Leverage = quantile(Leverage, 0.05),
            P95_Leverage = quantile(Leverage, 0.95), SD_Leverage = sd(Leverage),
            Median_AbsLeverage = median(AbsLeverage), P05_AbsLeverage = quantile(AbsLeverage, 0.05),
            P95_AbsLeverage = quantile(AbsLeverage, 0.95), SD_AbsLeverage = sd(AbsLeverage),
            Ref_Leverage_Real = first(Ref_Leverage_Real), Ref_Leverage_MC = first(Ref_Leverage_MC),
            Ref_AbsLeverage_Real = first(Ref_AbsLeverage_Real), Ref_AbsLeverage_MC = first(Ref_AbsLeverage_MC),
            .groups = "drop") %>%
  mutate(Pct_Bias_vs_MC     = (Median_Leverage - Ref_Leverage_MC) / abs(Ref_Leverage_MC) * 100,
         Abs_Bias_vs_MC     = Median_Leverage - Ref_Leverage_MC,
         Pct_Bias_vs_Real   = (Median_Leverage - Ref_Leverage_Real) / abs(Ref_Leverage_Real) * 100,
         Pct_Gap_Real_vs_MC = (Ref_Leverage_Real - Ref_Leverage_MC) / abs(Ref_Leverage_MC) * 100,
         AbsPct_Bias_vs_MC  = (Median_AbsLeverage - Ref_AbsLeverage_MC) / abs(Ref_AbsLeverage_MC) * 100,
         AbsAbs_Bias_vs_MC  = Median_AbsLeverage - Ref_AbsLeverage_MC)
write_csv(spatial_mc_summary, file.path(data_out_dir, "leverage_MC_summary_br.csv"))

bias_annot <- spatial_mc_summary %>% filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.2f pct-pts (%+.1f%%) bias", N, Abs_Bias_vs_MC, Pct_Bias_vs_MC))

p_leverage <- ggplot(spatial_mc_summary, aes(x = N, y = Median_Leverage)) +
  geom_ribbon(aes(ymin = P05_Leverage, ymax = P95_Leverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_Leverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_Leverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites, outlet included)", y = "Leverage (%, median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity (leverage) -- br",
       subtitle = "Dashed = reference leverage; dotted red = actual sample size") +
  theme_bw() +
  theme(legend.position = "bottom")
p_leverage
ggsave(file.path(plot_dir, "leverage_MC_plot_br.png"), p_leverage, width = 7, height = 5, dpi = 150)

bias_annot_abs <- spatial_mc_summary %>% filter(N == actual_n_sites) %>%
  mutate(label = sprintf("At N=%d: %+.2f pct-pts (%+.1f%%) bias", N, AbsAbs_Bias_vs_MC, AbsPct_Bias_vs_MC))

p_abs_leverage <- ggplot(spatial_mc_summary, aes(x = N, y = Median_AbsLeverage)) +
  geom_ribbon(aes(ymin = P05_AbsLeverage, ymax = P95_AbsLeverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_AbsLeverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_AbsLeverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  geom_vline(xintercept = actual_n_sites, color = "red", linetype = "dotted") +
  geom_text(data = bias_annot_abs, aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, size = 3.5, color = "red", inherit.aes = FALSE) +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Sample size (N sites, outlet included)", y = "|Leverage| (%, median, 5th-95th pct band)",
       title = "Spatial sample-size sensitivity (|leverage|) -- br",
       subtitle = "Dashed = reference |leverage|; dotted red = actual sample size") +
  theme_bw() +
  theme(legend.position = "bottom")
p_abs_leverage
ggsave(file.path(plot_dir, "abs_leverage_MC_plot_br.png"), p_abs_leverage, width = 7, height = 5, dpi = 150)


# compute optimums
bias_tol_abs   <- 5    # allowable |bias| in leverage percentage-points (vs. model truth) -- tune to taste
elbow_fraction <- 0.05

optimal_n_list <- list()
for (cc in constituents) {
  cc_summary <- spatial_mc_summary %>%
    filter(Constituent == cc) %>%
    arrange(N)
  
  n_vals        <- cc_summary$N
  width_pct     <- (cc_summary$P95_Leverage - cc_summary$P05_Leverage) / abs(cc_summary$Ref_Leverage_MC) * 100
  width_pct_abs <- (cc_summary$P95_AbsLeverage - cc_summary$P05_AbsLeverage) / abs(cc_summary$Ref_AbsLeverage_MC) * 100
  
  marginal_rate   <- -diff(width_pct) / diff(n_vals)
  elbow_idx       <- which(marginal_rate <= elbow_fraction * max(marginal_rate))[1]
  n_opt_precision <- if (is.na(elbow_idx)) NA_real_ else n_vals[elbow_idx + 1]
  
  marginal_rate_abs   <- -diff(width_pct_abs) / diff(n_vals)
  elbow_idx_abs       <- which(marginal_rate_abs <= elbow_fraction * max(marginal_rate_abs))[1]
  n_opt_precision_abs <- if (is.na(elbow_idx_abs)) NA_real_ else n_vals[elbow_idx_abs + 1]
  
  within_tol   <- abs(cc_summary$Abs_Bias_vs_MC) <= bias_tol_abs
  stays_within <- sapply(seq_along(within_tol), function(i) all(within_tol[i:length(within_tol)]))
  bias_idx     <- which(stays_within)[1]
  n_opt_bias   <- if (is.na(bias_idx)) NA_real_ else n_vals[bias_idx]
  
  within_tol_abs   <- abs(cc_summary$AbsAbs_Bias_vs_MC) <= bias_tol_abs
  stays_within_abs <- sapply(seq_along(within_tol_abs), function(i) all(within_tol_abs[i:length(within_tol_abs)]))
  bias_idx_abs     <- which(stays_within_abs)[1]
  n_opt_bias_abs   <- if (is.na(bias_idx_abs)) NA_real_ else n_vals[bias_idx_abs]
  
  optimal_n_list[[cc]] <- tibble(Constituent = cc,
                                 N_optimal_precision     = n_opt_precision,
                                 N_optimal_precision_abs = n_opt_precision_abs,
                                 N_optimal_bias          = n_opt_bias,
                                 N_optimal_bias_abs      = n_opt_bias_abs)
}
optimal_n <- bind_rows(optimal_n_list)
print(optimal_n)


#### Testing missing-data tolerance ####

missing_grid <- c(0.05, 0.10, 0.2, 0.3, 0.4, 0.50, 0.6, 0.7, 0.8)

missing_mc_results_list <- list()
for (cc in constituents) {
  cc_data <- clean_prepped %>%
    filter(!is.na(.data[[cc]])) %>%
    select(Site, CampaignID, Conc = all_of(cc), Area.m2, SpecificDischarge)
  
  campaign_outlet <- split(cc_data %>% filter(Site == outlet_site_id),
                           (cc_data %>% filter(Site == outlet_site_id))$CampaignID)
  campaign_sites  <- split(cc_data %>% filter(Site != outlet_site_id),
                           (cc_data %>% filter(Site != outlet_site_id))$CampaignID)
  
  ref_campaign_lev <- sapply(names(campaign_sites), function(camp) {
    sites  <- campaign_sites[[camp]]
    outlet <- campaign_outlet[[camp]]
    lev <- 100 * (sites$Conc - outlet$Conc) * (sites$Area.m2 / outlet$Area.m2) *
      (sites$SpecificDischarge / outlet$SpecificDischarge) / outlet$Conc
    c(Leverage = mean(lev), AbsLeverage = mean(abs(lev)))
  })
  ref_leverage_missing     <- mean(ref_campaign_lev["Leverage", ])
  ref_abs_leverage_missing <- mean(ref_campaign_lev["AbsLeverage", ])
  
  mc_results <- list()
  for (pct_missing in missing_grid) {
    lev_draws     <- numeric(n_iter)
    abs_lev_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      campaign_lev <- sapply(names(campaign_sites), function(camp) {
        sites  <- campaign_sites[[camp]]
        outlet <- campaign_outlet[[camp]]
        n_keep <- max(1, round(nrow(sites) * (1 - pct_missing)))
        kept   <- sites[sample(nrow(sites), n_keep), ]
        lev <- 100 * (kept$Conc - outlet$Conc) * (kept$Area.m2 / outlet$Area.m2) *
          (kept$SpecificDischarge / outlet$SpecificDischarge) / outlet$Conc
        c(Leverage = mean(lev), AbsLeverage = mean(abs(lev)))
      })
      lev_draws[i]     <- mean(campaign_lev["Leverage", ])
      abs_lev_draws[i] <- mean(campaign_lev["AbsLeverage", ])
    }
    mc_results[[paste(pct_missing)]] <- tibble(Constituent = cc, Pct_Missing = pct_missing * 100,
                                               Leverage = lev_draws, AbsLeverage = abs_lev_draws,
                                               Ref_Leverage_Missing = ref_leverage_missing,
                                               Ref_AbsLeverage_Missing = ref_abs_leverage_missing)
  }
  missing_mc_results_list[[cc]] <- bind_rows(mc_results)
}
missing_mc_results <- bind_rows(missing_mc_results_list)

missing_mc_summary <- missing_mc_results %>%
  group_by(Constituent, Pct_Missing) %>%
  summarize(Median_Leverage = median(Leverage), P05_Leverage = quantile(Leverage, 0.05),
            P95_Leverage = quantile(Leverage, 0.95), SD_Leverage = sd(Leverage),
            Median_AbsLeverage = median(AbsLeverage), P05_AbsLeverage = quantile(AbsLeverage, 0.05),
            P95_AbsLeverage = quantile(AbsLeverage, 0.95), SD_AbsLeverage = sd(AbsLeverage),
            Ref_Leverage_Missing = first(Ref_Leverage_Missing),
            Ref_AbsLeverage_Missing = first(Ref_AbsLeverage_Missing), .groups = "drop") %>%
  mutate(Pct_Bias    = (Median_Leverage - Ref_Leverage_Missing) / abs(Ref_Leverage_Missing) * 100,
         Abs_Bias    = Median_Leverage - Ref_Leverage_Missing,
         AbsPct_Bias = (Median_AbsLeverage - Ref_AbsLeverage_Missing) / abs(Ref_AbsLeverage_Missing) * 100,
         AbsAbs_Bias = Median_AbsLeverage - Ref_AbsLeverage_Missing)
write_csv(missing_mc_summary, file.path(data_out_dir, "leverage_missing_data_summary_br.csv"))

p_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_Leverage)) +
  geom_ribbon(aes(ymin = P05_Leverage, ymax = P95_Leverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_Leverage_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "Leverage (%, median, 5th-95th pct band)",
       title = "Missing-data tolerance (leverage) -- br",
       subtitle = "Blue dashed = no-missing-data baseline (3 real campaigns); outlet never dropped") +
  theme_bw()
p_missing
ggsave(file.path(plot_dir, "leverage_missing_data_plot_br.png"), p_missing, width = 7, height = 5, dpi = 150)

p_abs_missing <- ggplot(missing_mc_summary, aes(x = Pct_Missing, y = Median_AbsLeverage)) +
  geom_ribbon(aes(ymin = P05_AbsLeverage, ymax = P95_AbsLeverage), fill = "grey70", alpha = 0.5) +
  geom_line() + geom_point(size = 1.5) +
  geom_hline(aes(yintercept = Ref_AbsLeverage_Missing), color = "blue", linetype = "dashed") +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = "Missing data (%)", y = "|Leverage| (%, median, 5th-95th pct band)",
       title = "Missing-data tolerance (|leverage|) -- br",
       subtitle = "Blue dashed = no-missing-data baseline (3 real campaigns); outlet never dropped") +
  theme_bw()
p_abs_missing
ggsave(file.path(plot_dir, "abs_leverage_missing_data_plot_br.png"), p_abs_missing, width = 7, height = 5, dpi = 150)

### calculate optimums 

bias_tol_missing_abs  <- 1      # leverage percentage-point bias tolerance -- tune to taste
width_tol_missing_abs <- 5      # leverage percentage-point width tolerance -- tune to taste

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
    mutate(Width_abs     = P95_Leverage - P05_Leverage,
           Width_abs_Lev = P95_AbsLeverage - P05_AbsLeverage)
  
  pct_allowable_bias         <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$Abs_Bias), bias_tol_missing_abs)
  pct_allowable_width        <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_abs, width_tol_missing_abs)
  pct_allowable_bias_absLev  <- find_allowable_pct(cc_summary$Pct_Missing, abs(cc_summary$AbsAbs_Bias), bias_tol_missing_abs)
  pct_allowable_width_absLev <- find_allowable_pct(cc_summary$Pct_Missing, cc_summary$Width_abs_Lev, width_tol_missing_abs)
  
  allowable_missing_list[[cc]] <- tibble(Constituent = cc,
                                         Allowable_Pct_Missing_Bias         = pct_allowable_bias,
                                         Allowable_Pct_Missing_Width        = pct_allowable_width,
                                         #Allowable_Pct_Missing              = min(pct_allowable_bias, pct_allowable_width),
                                         Allowable_Pct_Missing_Bias_AbsLev  = pct_allowable_bias_absLev,
                                         Allowable_Pct_Missing_Width_AbsLev = pct_allowable_width_absLev)
  #Allowable_Pct_Missing_AbsLev       = min(pct_allowable_bias_absLev, pct_allowable_width_absLev))
}
allowable_missing <- bind_rows(allowable_missing_list)
print(allowable_missing)

#
#### Varying spatial distribution of sites (mainstem vs. tributary) ####

shift_fraction <- 0.99

total_watershed_area  <- max(predpts_attrs$Area_m2)
mainstem_frac_of_total <- 0.06
mainstem_threshold     <- mainstem_frac_of_total * total_watershed_area

real_is_mainstem <- area_lookup[real_predids] >= mainstem_threshold
p_real <- mean(real_is_mainstem)

nn_dist_mat <- as.matrix(dist(predpts_attrs[, c("X", "Y")]))
diag(nn_dist_mat) <- NA
print(summary(apply(nn_dist_mat, 1, min, na.rm = TRUE)))

min_separation_mainstem_m <- 350
min_separation_trib_m     <- 400

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
mainstem_predids_pool <- setdiff(names(area_lookup)[is_mainstem_all], outlet_predID)
trib_predids_pool     <- setdiff(names(area_lookup)[!is_mainstem_all], outlet_predID)

mainstem_predids_thin <- decluster_points(mainstem_predids_pool, predpts_attrs, min_separation_mainstem_m, seed = 1)
trib_predids_thin     <- decluster_points(trib_predids_pool,     predpts_attrs, min_separation_trib_m, seed = 2)

cat("Mainstem pool:", length(mainstem_predids_pool), "-> thinned to", length(mainstem_predids_thin), "\n")
cat("Tributary pool:", length(trib_predids_pool), "-> thinned to", length(trib_predids_thin), "\n")
cat("Need at least", actual_n_sites - 1, "non-outlet candidates available in the more restrictive pool.\n")

predpts_attrs <- predpts_attrs %>%
  mutate(
    Group = case_when(
      predID == outlet_predID           ~ "Outlet",
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
    "Outlet"                     = "black",
    "Mainstem (sampling pool)"  = "red",
    "Tributary (sampling pool)" = "orange",
    "Mainstem (thinned out)"    = "grey",
    "Tributary (thinned out)"   = "#2e86c1"
  )) +
  labs(title = paste0(toupper("br"), ": mainstem/tributary sampling pools"),
       x = NULL, y = NULL, color = NULL) +
  theme_minimal()
p_network_map

distribution_mc_results_list <- list()
for (cc in constituents) {
  fs      <- field_setups[[cc]]
  m       <- fs$setup$m
  mu_cond <- fs$setup$mu_cond
  L_space <- fs$setup$L_space
  
  area_by_pred <- area_lookup[paste0("P", seq_len(m))]
  is_mainstem  <- area_by_pred >= mainstem_threshold
  outlet_idx   <- which(paste0("P", seq_len(m)) == outlet_predID)
  mainstem_idx <- match(mainstem_predids_thin, paste0("P", seq_len(m)))
  trib_idx     <- match(trib_predids_thin,     paste0("P", seq_len(m)))
  
  p_real <- mean(is_mainstem[real_predids])
  
  scenarios <- tibble(
    Scenario = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"),
    Target_mainstem_frac = c(p_real * (1 - shift_fraction), p_real,
                             p_real + shift_fraction * (1 - p_real))
  )
  
  z             <- matrix(rnorm(m * n_iter), nrow = m, ncol = n_iter)
  spatial_noise <- t(L_space) %*% z
  sim_mat       <- exp(mu_cond + spatial_noise)
  
  ref_leverage_real     <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_Leverage_Real)
  ref_abs_leverage_real <- ref_leverage_real_by_cc %>% filter(Constituent == cc) %>% pull(Ref_AbsLeverage_Real)
  
  full_pop_leverage <- apply(sim_mat, 2, function(vals) {
    lev_pct <- 100 * (vals - vals[outlet_idx]) * (area_by_pred / area_by_pred[outlet_idx]) / vals[outlet_idx]
    c(Leverage = mean(lev_pct[-outlet_idx]), AbsLeverage = mean(abs(lev_pct[-outlet_idx])))
  })
  ref_leverage_mc     <- mean(full_pop_leverage["Leverage", ])
  ref_abs_leverage_mc <- mean(full_pop_leverage["AbsLeverage", ])
  
  mc_results <- list()
  for (s in seq_len(nrow(scenarios))) {
    scen_name  <- scenarios$Scenario[s]
    n_mainstem <- min(max(round(scenarios$Target_mainstem_frac[s] * (actual_n_sites - 1)), 0), actual_n_sites - 1)
    n_trib     <- (actual_n_sites - 1) - n_mainstem
    
    lev_draws     <- numeric(n_iter)
    abs_lev_draws <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      idx     <- c(outlet_idx, sample(mainstem_idx, n_mainstem), sample(trib_idx, n_trib))
      vals    <- sim_mat[idx, i]
      areas   <- area_by_pred[idx]
      lev_pct <- 100 * (vals - vals[1]) * (areas / areas[1]) / vals[1]
      lev_draws[i]     <- mean(lev_pct[-1])
      abs_lev_draws[i] <- mean(abs(lev_pct[-1]))
    }
    mc_results[[scen_name]] <- tibble(Constituent = cc, Scenario = scen_name,
                                      Leverage = lev_draws, AbsLeverage = abs_lev_draws,
                                      Ref_Leverage_Real = ref_leverage_real, Ref_Leverage_MC = ref_leverage_mc,
                                      Ref_AbsLeverage_Real = ref_abs_leverage_real, Ref_AbsLeverage_MC = ref_abs_leverage_mc)
  }
  distribution_mc_results_list[[cc]] <- bind_rows(mc_results)
}
distribution_mc_results <- bind_rows(distribution_mc_results_list)

distribution_mc_summary <- distribution_mc_results %>%
  mutate(Scenario = factor(Scenario, levels = c("Tributary-dominated", "Uniform (real mix)", "Mainstem-dominated"))) %>%
  group_by(Constituent, Scenario) %>%
  summarize(Median_Leverage = median(Leverage), P05_Leverage = quantile(Leverage, 0.05),
            P95_Leverage = quantile(Leverage, 0.95), SD_Leverage = sd(Leverage),
            Median_AbsLeverage = median(AbsLeverage), P05_AbsLeverage = quantile(AbsLeverage, 0.05),
            P95_AbsLeverage = quantile(AbsLeverage, 0.95), SD_AbsLeverage = sd(AbsLeverage),
            Ref_Leverage_Real = first(Ref_Leverage_Real), Ref_Leverage_MC = first(Ref_Leverage_MC),
            Ref_AbsLeverage_Real = first(Ref_AbsLeverage_Real), Ref_AbsLeverage_MC = first(Ref_AbsLeverage_MC),
            .groups = "drop") %>%
  mutate(Abs_Bias_vs_MC        = Median_Leverage - Ref_Leverage_MC,
         Abs_Bias_vs_Real      = Median_Leverage - Ref_Leverage_Real,
         AbsAbs_Bias_vs_MC     = Median_AbsLeverage - Ref_AbsLeverage_MC,
         AbsAbs_Bias_vs_Real   = Median_AbsLeverage - Ref_AbsLeverage_Real,
         Width_abs             = P95_Leverage - P05_Leverage,
         Width_abs_AbsLeverage = P95_AbsLeverage - P05_AbsLeverage)
write_csv(distribution_mc_summary, file.path(data_out_dir, "leverage_distribution_summary_br.csv"))

p_dist_leverage <- ggplot(distribution_mc_summary, aes(x = Scenario, y = Median_Leverage)) +
  geom_pointrange(aes(ymin = P05_Leverage, ymax = P95_Leverage)) +
  geom_hline(aes(yintercept = Ref_Leverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_Leverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Leverage (%, median, 5th-95th pct band)",
       title = "Spatial sampling-distribution sensitivity (leverage) -- br",
       subtitle = paste0("N fixed at actual_n_sites = ", actual_n_sites, " (outlet included); dashed = reference leverage")) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
p_dist_leverage
ggsave(file.path(plot_dir, "leverage_distribution_plot_br.png"), p_dist_leverage, width = 7, height = 5, dpi = 150)

p_abs_dist_leverage <- ggplot(distribution_mc_summary, aes(x = Scenario, y = Median_AbsLeverage)) +
  geom_pointrange(aes(ymin = P05_AbsLeverage, ymax = P95_AbsLeverage)) +
  geom_hline(aes(yintercept = Ref_AbsLeverage_MC, color = "Model truth (MC population)"), linetype = "dashed") +
  geom_hline(aes(yintercept = Ref_AbsLeverage_Real, color = "Realized truth (3 real campaigns)"), linetype = "dashed") +
  scale_color_manual(name = "Reference",
                     values = c("Model truth (MC population)" = "blue",
                                "Realized truth (3 real campaigns)" = "darkgreen")) +
  facet_wrap(~ Constituent, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "|Leverage| (%, median, 5th-95th pct band)",
       title = "Spatial sampling-distribution sensitivity (|leverage|) -- br",
       subtitle = paste0("N fixed at actual_n_sites = ", actual_n_sites, " (outlet included); dashed = reference |leverage|")) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))
p_abs_dist_leverage
ggsave(file.path(plot_dir, "abs_leverage_distribution_plot_br.png"), p_abs_dist_leverage, width = 7, height = 5, dpi = 150)

#