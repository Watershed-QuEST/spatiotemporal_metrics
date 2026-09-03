#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-08-23 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-02

# Requires: 
# 1. 02_build_synthetic_data.R must be run first to produce data/nm_clean.csv, data/nm_field_setups.rds, and data/nm_synthetic_extended.csv. If these are available, no need to rerun 02_build_synthetic_data.R
# 2. spatiotemporal_helpers.R must be in the same folder as this script.

# This script demonstrates how to calculate CVt, the temporal coefficient of variation (across sampling campaigns, one value per sampling site). 
# It also performs a sensitivity analysis to answer the question: how much does a CVt estimate vary/undershoot depending how many sampling campaigns are done? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possible fits to a SSN2 spatial covariance model. 
# Ref_CVt (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R, which simulates 60 (1 per month x 5 years) synoptic campaigns.

# This script is organized as:
#   PART A -- CVt_observed: computed directly from the real nm data
#   PART B -- Monte Carlo sensitivity sweep: 

# To adapt this script to a metric other than CVt: swap out the statistic computed inside the two n_iter loops in Part B, and Part A's real-data calculation to match.

# Outputs:
#   data/CVt_observed.csv: real CVt per constituent
#   data/CVt_MC_summary.csv: sensitivity analysis summary (median, P05-P95, bias vs Ref_CVt)
#   plots/CVt_MC_plot.png

#### Packages ####
library(tidyverse)
source("spatiotemporal_helpers.R")

#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Monte Carlo iterations
n_iter    <- 10000                                    
camp_grid <- c(2, 3, 4, 5, 7, 10, 15, 20, 30, 50)      # campaign counts to test

n_synthetic_months <- 60      # must match 02_build_synthetic_data.R
time_unit_days     <- 30.44   # must match 02_build_synthetic_data.R
total_days <- n_synthetic_months * time_unit_days  # span campaign times are drawn from

set.seed(42)

clean = read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups = readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended = read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)

#### PART A -- Calculate CVt of real toy dataset ####
# Per-site CV across that site's real campaigns

# for each site and constituent
CVt_observed <- map_dfr(names(field_setups), function(cc) {
  per_site <- clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(Site) %>%
    summarize(n_campaigns = n(), CVt = sd(.data[[cc]]) / mean(.data[[cc]]), .groups = "drop") %>%
    filter(n_campaigns >= 3)
  tibble(Constituent = cc, n_sites_used = nrow(per_site), CVt_observed = per_site$CVt)
})
print(CVt_observed)

# averaged across sites for each constituent
CVt_observed_mean <- map_dfr(names(field_setups), function(cc) {
  per_site <- clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(Site) %>%
    summarize(n_campaigns = n(), CVt = sd(.data[[cc]]) / mean(.data[[cc]]), .groups = "drop") %>%
    filter(n_campaigns >= 3)
  tibble(Constituent = cc, n_sites_used = nrow(per_site), CVt_observed = mean(per_site$CVt))
})
print(CVt_observed_mean)
