#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-07

# Requires: 02_build_synthetic_data.R must be run to produce data/[dataset]_clean.csv for each dataset. 

# This script demonstrates how to calculate CVs, the spatial coefficient of variation (across sites, one value per sampling event) on three datasets: NM, BR, and Gu et al. 2021. 
# Complementary sensitivity analyses can be found in CVs_sensitivity.R

# Outputs:
#   data/CVs_by_campaign.csv:  real CVs, one row per real campaign x constituent
#   data/CVs_observed.csv: Mean_CVs and SD_CVs across those real campaigns

#### Packages ####
library(tidyverse)

#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

clean            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)

#### Calculate CVs of real toy dataset ####
# Each site's concentration averaged across its available real campaigns, then CV taken across sites.

## For each constituent in `constituents`: CVs computed separately for each real campaign (spatial CV across whichever sites were sampled that campaign). Campaigns with fewer than 2 non-NA sites get CVs = NA (sd() of a single value is undefined).

CVs_by_campaign_list <- list()
CVs_observed_list    <- list()

for (cc in constituents) {
  
  by_campaign <- clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(CampaignID) %>%
    summarize(n_sites = n(),
              CVs = if (n() >= 2) sd(.data[[cc]]) / mean(.data[[cc]]) else NA_real_,
              .groups = "drop") %>%
    mutate(Constituent = cc, .before = 1)
  
  cat(cc, "-- observed CVs per real campaign:\n")
  print(by_campaign)
  
  ## Summarized across campaigns - Mean_CVs is their average
  valid_cvs <- by_campaign$CVs[!is.na(by_campaign$CVs)]
  observed  <- tibble(Constituent = cc, n_campaigns = length(valid_cvs),
                      Mean_CVs = mean(valid_cvs),
                      SD_CVs = if (length(valid_cvs) >= 2) sd(valid_cvs) else NA_real_)
  
  CVs_by_campaign_list[[cc]] <- by_campaign
  CVs_observed_list[[cc]]    <- observed
}

CVs_by_campaign <- bind_rows(CVs_by_campaign_list)
write_csv(CVs_by_campaign, file.path(data_out_dir, "CVs_by_campaign.csv"))

CVs_observed <- bind_rows(CVs_observed_list)
write_csv(CVs_observed, file.path(data_out_dir, "CVs_observed.csv"))