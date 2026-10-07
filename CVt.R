#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-08-23 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-08
# Bre Rivera Waterman, 2026-09-09 pulling in package/helper function and comparing to previous calculations
# Bre Rivera Waterman,09/29/2026 for conceptual fig values
# Alex Webster, 2026-10-06 to add br analysis and remove sensitivity analysis


# Requires: 02_build_synthetic_data.R must be run first to produce data/nm_clean.csv, data/nm_field_setups.rds, and data/nm_synthetic_extended.csv. If these are available, no need to rerun 02_build_synthetic_data.R

# This script demonstrates how to calculate CVt, the temporal coefficient of variation (across sampling campaigns, one value per sampling site). 
# It also performs a sensitivity analysis to answer the question: how much does a CVt estimate vary/undershoot depending how many sampling campaigns are done? This analysis uses Monte Carlo iteration so the resulting uncertainty bands reflect fresh draws from all possible fits to a SSN2 spatial covariance model. 
# Ref_CVt (the "population truth" to compare against) comes from the extended static dataset created in 02_build_synthetic_data.R, which simulates 60 (1 per month x 5 years) synoptic campaigns.

# This script demonstrates how to calculate CVs, the spatial coefficient of variation (across sites, one value per sampling event) on three datasets: NM, BR, and Gu et al. 2021. 
# Complementary sensitivity analyses can be found in CVs_sensitivity.R

# Outputs:
#   data/CVt_observed.csv:      real CVt per catchment x site x constituent
#   data/CVt_observed_mean.csv: real CVt averaged across catchment x sites, per constituent

#### Packages ####
library(tidyverse)
source("cv_fxhelper.R")


#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

clean_nm            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
clean_nm$Catchment = "nm"
clean_br            <- read_csv(file.path(data_out_dir, "br_clean.csv"), show_col_types = FALSE)
clean_br$Catchment = "br"

clean <- bind_rows(clean_nm, clean_br)

actual_n_campaigns <- n_distinct(clean$CampaignID)  # real number of campaigns conducted across catchments
actual_n_campaigns_nm <- n_distinct(clean_nm$CampaignID)  # real number of nm campaigns conducted
actual_n_campaigns_br <- n_distinct(clean_br$CampaignID)  # real number of br campaigns conducted


#### PART A.1 -- Calculate CVt of real toy dataset ####
# Per-site CV across that site's real campaigns

## For each constituent in `constituents`: per-site CV across that site's real campaigns (sites with fewer than 3 real campaigns are dropped, since sd() of 1-2 points is too noisy to call a temporal CV).
CVt_observed_list      <- list()
CVt_observed_mean_list <- list()

for (cc in constituents) {
  
  per_site_cvt <- clean %>%
    filter(!is.na(.data[[cc]])) %>%
    group_by(Catchment, Site) %>%
    summarize(n_sites = n(), CVt = sd(.data[[cc]]) / mean(.data[[cc]]), .groups = "drop") %>%
    filter(n_sites >= 3)

  print(per_site_cvt)
  
  CVt_observed_list[[cc]] <- per_site_cvt %>%
    mutate(Constituent = cc, .before = 1)
  
  CVt_observed_mean_list[[cc]] <- per_site_cvt %>%
    group_by(Catchment) %>%
    summarize(n_sites_used = n(),
              CVt_observed = mean(CVt),
              .groups = "drop") %>%
    mutate(Constituent = cc, .before = 1)
}


## ---- Combine and save ----
CVt_observed <- bind_rows(CVt_observed_list)
write_csv(CVt_observed, file.path(data_out_dir, "CVt_observed.csv"))

CVt_observed_mean <- bind_rows(CVt_observed_mean_list)
write_csv(CVt_observed_mean, file.path(data_out_dir, "CVt_observed_mean.csv"))
print(CVt_observed_mean)


#### PART A.2 -- Calculate CVt of real toy dataset using package!! ####
# CVt is calculated for one site across campaigns.

CVt_result <- temporal_cv(data = clean,
                         concentration = constituents,
                         site = "Site",
                         watershed = "Catchment",
                         digits = 2)

# One row per site × constituent
CVt_by_site.2 <- 
  CVt_result$by_site %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            Site,
            n_sites = n_used,
            CVt = temporal_cv)

print(CVt_by_site.2)

write_csv(CVt_by_site.2, file.path(data_out_dir, "CVt_by_site.2.csv"))

# Mean and SD across campaign-level CVs for each constituent
CVt_observed.2 <- CVt_result$watershed_summary %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            n_sites = n_sites_used,
            Mean_CVt = mean_temporal_cv,
            SD_CVt = sd_temporal_cv )

print(CVt_observed.2)

write_csv(CVt_observed.2, file.path(data_out_dir, "CVt_observed.2.csv"))

#### PART A.3 -- Compare original and package results ####
site_comparison <- CVt_observed %>%
  mutate(CVt_original_rounded = round(CVt, 2)) %>%
  full_join(CVt_by_site.2 %>%
              rename(
                #n_sites_package = n_sites,
                     CVt_package = CVt ),
            by = c("Catchment","Constituent", "Site")) %>%
  #rename(n_sites_original = n_sites) %>%
  mutate(#n_sites_difference = n_sites_package - n_sites_original,
         CVt_difference = CVt_package - CVt_original_rounded) %>%
  arrange(Constituent, Site)

print(site_comparison)


#watershed summaries
original_by_catchment <- CVt_observed %>%
  group_by(Constituent, Catchment) %>%
  summarize(n_sites_original = n(),
            Mean_CVt_original = mean(CVt),
            SD_CVt_original = if (n() >= 2) sd(CVt) else NA_real_,
            .groups = "drop") %>%
  mutate(Mean_CVt_original_rounded = round(Mean_CVt_original, 2),
         SD_CVt_original_rounded = round(SD_CVt_original, 2))

package_by_catchment <- CVt_observed.2 %>%
  rename(n_sites_package = n_sites,
         Mean_CVt_package = Mean_CVt,
         SD_CVt_package = SD_CVt)

summary_comparison <- original_by_catchment %>%
  full_join(package_by_catchment, by = c("Constituent", "Catchment")) %>%
  mutate(Mean_CVt_difference = Mean_CVt_package - Mean_CVt_original_rounded,
         SD_CVt_difference = SD_CVt_package - SD_CVt_original_rounded) %>%
  arrange(Constituent, Catchment)

print(summary_comparison, n = Inf)


#visual comparison
CVt_plot_data <- 
  bind_rows(CVt_observed %>%
              transmute(Catchment, Constituent, Site, CVt,
                        Approach = "Original QuEST calculation" ),
            CVt_by_site.2 %>%
              transmute(Catchment, Constituent,Site, CVt, Approach = "CV helper")) %>%
  filter(!is.na(CVt)) %>%
  mutate(Approach = factor(Approach,
                           levels = c("Original QuEST calculation", "CV helper") ) )


p_cvt_comparison <- 
  ggplot(CVt_plot_data, aes(x = Approach, y = CVt, fill = Approach)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    color = NA ) +
  geom_jitter(
    aes(color = Catchment),
    width = 0.05,
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
  labs(title = "Temporal CV by calculation approach",
       subtitle = "Each point is one sampling location",
       x = NULL,
       y = "Temporal coefficient of variation (CVt)",
       fill = "Approach") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

p_cvt_comparison


#### PART A.4 -- overall average for conceptual figure ####
CVt_plot_data %>% group_by(Catchment, Constituent) %>% summarise(mean = mean(CVt), sd = sd(CVt))
