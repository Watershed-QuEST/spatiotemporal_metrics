#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-07-28 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-07
# Bre Rivera Waterman, 2026-09-01 pulling in package/helper function and comparing to previous calculations

# Requires: 02_build_synthetic_data.R must be run to produce data/[dataset]_clean.csv for each dataset. 

# This script demonstrates how to calculate CVs, the spatial coefficient of variation (across sites, one value per sampling event) on three datasets: NM, BR, and Gu et al. 2021. 
# Complementary sensitivity analyses can be found in CVs_sensitivity.R

# Outputs:
#   data/CVs_by_campaign.csv:  real CVs, one row per real campaign x constituent
#   data/CVs_observed.csv: Mean_CVs and SD_CVs across those real campaigns

#### Packages ####
library(tidyverse)
source("cv_helper.R")

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

#### PART A.1 -- Calculate CVs of real toy dataset ####
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


#### PART A.2 -- Calculate CVs of real toy dataset using package!! ####
# CVs is calculated across sites within each campaign.

CVs_result <- spatial_cv(data = clean,
                         concentration = names(field_setups),
                         event = "CampaignID",
                         digits = 2)

# One row per campaign × constituent
CVs_by_campaign.2 <- 
  CVs_result$by_event %>%
  transmute(Constituent = constituent,
            CampaignID,
            n_sites = n_used,
            CVs = spatial_cv)

print(CVs_by_campaign.2)

write_csv(CVs_by_campaign.2, file.path(data_out_dir, "CVs_by_campaign.2.csv"))

# Mean and SD across campaign-level CVs for each constituent
CVs_observed.2 <- CVs_result$watershed_summary %>%
  transmute(Constituent = constituent,
            n_campaigns = n_events_used,
            Mean_CVs = mean_spatial_cv,
            SD_CVs = sd_spatial_cv )

print(CVs_observed.2)

write_csv(CVs_observed.2, file.path(data_out_dir, "CVs_observed.2.csv"))



#### PART A.3 -- Compare original and package results ####

campaign_comparison <- CVs_by_campaign %>%
  mutate(CVs_original_rounded = round(CVs, 2)) %>%
  full_join(CVs_by_campaign.2 %>%
              rename(n_sites_package = n_sites,
                     CVs_package = CVs ),
            by = c("Constituent", "CampaignID")) %>%
  rename(n_sites_original = n_sites) %>%
  mutate(n_sites_difference = n_sites_package - n_sites_original,
         CVs_difference = CVs_package - CVs_original_rounded) %>%
  arrange(Constituent, CampaignID)

print(campaign_comparison, n = Inf)


#watershed summaries
summary_comparison <- CVs_observed %>%
  mutate(Mean_CVs_original_rounded = round(Mean_CVs, 2),
         SD_CVs_original_rounded = round(SD_CVs, 2)) %>%
  full_join(CVs_observed.2 %>%
              rename(
                n_campaigns_package = n_campaigns,
                Mean_CVs_package = Mean_CVs,
                SD_CVs_package = SD_CVs),
            by = "Constituent") %>%
  rename(n_campaigns_original = n_campaigns) %>%
  mutate(n_campaigns_difference = n_campaigns_package - n_campaigns_original,
         Mean_CVs_difference = Mean_CVs_package - Mean_CVs_original_rounded,
         SD_CVs_difference = SD_CVs_package - SD_CVs_original_rounded) %>%
  arrange(Constituent)

print(summary_comparison, n = Inf)

#visual comparison
CVs_plot_data <- 
  bind_rows(CVs_by_campaign %>%
              transmute(Constituent, CampaignID, CVs,
                        Approach = "Original QuEST calculation" ),
            CVs_by_campaign.2 %>%
              transmute(Constituent,CampaignID, CVs, Approach = "CV helper")) %>%
  filter(!is.na(CVs)) %>%
  mutate(Approach = factor(Approach,
                           levels = c("Original QuEST calculation", "CV helper") ) )


p_cv_comparison <- 
  ggplot(CVs_plot_data, aes(x = Approach, y = CVs, fill = Approach)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    color = NA ) +
  geom_jitter(
    width = 0.07,
    height = 0,
    size = 1.8,
    alpha = 0.75 ) +
  facet_wrap(~ Constituent,
             scales = "free_y",
             ncol = 2) +
  scale_fill_manual(
    values = c(
      "Original QuEST calculation" = "#4C78A8",
      "CV helper" = "#F58518" )) +
  labs(title = "Spatial CV by calculation approach",
       subtitle = "Each point is one sampling campaign",
       x = NULL,
       y = "Spatial coefficient of variation (CVs)",
       fill = "Approach") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

p_cv_comparison
