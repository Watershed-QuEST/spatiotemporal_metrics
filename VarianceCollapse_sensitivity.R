#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Joanna Blaszczak, 2026-11-06 (adapted from Webster CVs_sensitivity)
# Last update (Person, Date): Joanna Blaszczak, 2026-11-06

# Requires: 02_build_synthetic_data.R must be run first (produces data/*_clean.csv, data/*_field_setups.rds, and data/*_synthetic_extended.csv).

# This script performs a sensitivity analysis to answer the questions 
# 1) how much does a variance collapse (VC) estimate vary/undershoot
# depending on how many sites you sample? This analysis uses Monte Carlo
# iteration so the resulting uncertainty bands reflect fresh draws from all possibile fits to a SSN2 spatial covariance model. 
# Ref_CVs (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R,
# which simulates a synoptic campaign that sampled every 250 m along the network.

# Outputs:
#   data/VarianceCollapse_MC_summary.csv:  MC summary
#   plots/VarianceCollapse_MC_plot.png: VC vs N, with reference line and bias annotation
#   plots/VarianceCollapse_sensitivity_metrics_plot.png: percent bias and SD vs N, side by side



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













