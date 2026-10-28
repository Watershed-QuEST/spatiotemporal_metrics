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
