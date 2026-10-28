#### READ ME ####

# Project: QuEST Spatiotemporal Metrics Commentary
# Author: Alex Webster, 2026-08-23 (with help building complex helper functions from Claude version 1.24012.9 (03c61d) 2026-07-24T04:59:17.000Z... heavily reviewed and edited by A. Webster)
# Last update (Person, Date): Alex Webster, 2026-09-08
# Bre Rivera Waterman, 2026-09-09 pulling in package/helper function and comparing to previous calculations
# Bre Rivera Waterman,09/29/2026 for conceptual fig values
# Alex Webster, 2026-10-06 to add br analysis and remove sensitivity analysis
# Alex Webster, 2026-10-08 to add extended synthetic data to make CVt n (3 per watershed) closer to CVs n (~20) for comparison


# Requires: 02_build_synthetic_data.R must be run first to produce data/[dataset]_clean.csv and data/[dataset]_synthetic_extended.csv. If these are available, no need to rerun 02_build_synthetic_data.R

# This script demonstrates how to calculate CVt, the temporal coefficient of variation (across sampling campaigns, one value per sampling site) on three datasets: NM, BR, and Gu et al. 2021. 
# Complementary sensitivity analyses can be found in CVt_sensitivity.R

# Outputs:
#   data/CVt_observed.csv:      real CVt per catchment x site x constituent
#   data/CVt_observed_mean.csv: real CVt averaged across catchment x sites, per constituent
#   data/CVt_observed_extended.csv:      CVt using extended synthetic data per catchment x site x constituent
#   data/CVt_observed_extended_mean.csv: CVt using extended synthetic data, averaged across catchment x sites, per constituent

#### Packages ####
library(tidyverse)
source("cv_fxhelper.R")
library(SSN2)

#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

# chem data from 3 consecutive real campaigns in nm and br
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


## ---- Combine and save
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
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p_cvt_comparison

#### PART A.4 -- calculate using extended synthetic data ####

### Load data
# chem data from synthetically extended data to use real sites but extend to 20 campaigns
# nm
extended_nm            <- read_csv(file.path(data_out_dir, "nm_synthetic_with_real.csv"), show_col_types = FALSE)
# the predpts that actually have a real-data equivalent -- these are "the real sites"
real_predIDs_nm <- extended_nm %>%
  filter(Source == "real") %>%
  distinct(predID) %>%
  pull(predID)
# 20 consecutive campaigns, beginning at the first campaign with real data
start_campaign_nm <- extended_nm %>%
  filter(Source == "real") %>%
  pull(CampaignNum) %>%
  min()
campaign_range_nm <- start_campaign_nm:(start_campaign_nm + 19)
real_sites_extended_nm <- extended_nm %>%
  filter(predID %in% real_predIDs_nm, CampaignNum %in% campaign_range_nm)
# br
extended_br            <- read_csv(file.path(data_out_dir, "br_synthetic_with_real.csv"), show_col_types = FALSE)
extended <- bind_rows(extended_nm, extended_br)
# the predpts that actually have a real-data equivalent -- these are "the real sites"
real_predIDs_br <- extended_br %>%
  filter(Source == "real") %>%
  distinct(predID) %>%
  pull(predID)
# 20 consecutive campaigns, beginning at the first campaign with real data
start_campaign_br <- extended_br %>%
  filter(Source == "real") %>%
  pull(CampaignNum) %>%
  min()
campaign_range_br <- start_campaign_br:(start_campaign_br + 19)
real_sites_extended_br <- extended_br %>%
  filter(predID %in% real_predIDs_br, CampaignNum %in% campaign_range_br)

# combine
real_sites_extended = bind_rows(real_sites_extended_br, real_sites_extended_nm)
# pivot
real_sites_extended_wide <- real_sites_extended %>%
  pivot_wider(
    names_from  = Constituent,
    values_from = c(Conc, Source)
  )
real_sites_extended_wide <- real_sites_extended_wide %>%
  rename_with(~ str_remove(., "\\.\\.mg\\..*$"))

## calculate CVt using package
CVt_extended_result <- temporal_cv(
  data = real_sites_extended_wide,
  watershed = "Catchment",
  concentration = c("Conc_NPOC", "Conc_TDN"),
  site = "predID",
  digits = 2)

# One row per site × constituent
CVt_extended_by_site.2 <- 
  CVt_extended_result$by_site %>%
  transmute(Catchment = Catchment,
            Constituent = constituent,
            predID,
            n_sites = n_used,
            CVt = temporal_cv)

print(CVt_extended_by_site.2)

write_csv(CVt_extended_by_site.2, file.path(data_out_dir, "CVt_extended_by_site.csv"))



#### PART A.5 -- overall average for conceptual figure ####

# real data only
CVt_plot_data %>% group_by(Catchment, Constituent) %>% summarise(mean = mean(CVt), sd = sd(CVt))

# real + synthetic data
CVt_extended_by_site.2 %>% group_by(Catchment, Constituent) %>% summarise(mean = mean(CVt), sd = sd(CVt))


#### PART B -- Final demonstration figures ####

## QuEST calculation compared to helper/package function
p_cvt_comparison
ggsave(file.path(plot_dir, "CVt_packagetoQuEST_comparision.png"), p_cvt_comparison, width = 9, height = 6, dpi = 150)

## CVs vs CVt: 

# load CVs data
CVs_by_campaign = read_csv(file.path(data_out_dir, "CVs_by_campaign.2.csv"))
CVs_extended_by_campaign = read_csv(file.path(data_out_dir, "CVs_extended_by_campaign.csv"))

# Max CV
cv_max <- function(n) sqrt(n)

cv_n <- tribble(
  ~Catchment, ~n_sites, ~n_campaigns,
  "nm",       22,       3,
  "br",       19,       3
)

cv_n_extended <- tribble(
  ~Catchment, ~n_sites, ~n_campaigns,
  "nm",       22,       20,
  "br",       19,       20
)

cv_plot_df <- bind_rows(
  biplot_real_df %>% mutate(Dataset = "Real only")                    %>% left_join(cv_n,          by = "Catchment"),
  biplot_ext_df  %>% mutate(Dataset = "Extended (real + synthetic)")  %>% left_join(cv_n_extended, by = "Catchment")
) %>%
  mutate(
    CVs_max       = cv_max(n_sites),
    CVt_max       = cv_max(n_campaigns),
    Mean_CVs_frac = mean_CVs / CVs_max,
    Mean_CVt_frac = mean_CVt / CVt_max,
    SD_CVs_frac   = sd_CVs / CVs_max,
    SD_CVt_frac   = sd_CVt / CVt_max
  )

#Combine with the mean CVt/CVs biplot data
axis_limit      <- 1.5
raw_slope       <- 1
max_ratio_slope <- mean(cv_plot_df$CVt_max) / mean(cv_plot_df$CVs_max)

max_ratio_by_dataset <- cv_plot_df %>%
  group_by(Dataset) %>%
  summarize(max_ratio_slope = mean(CVt_max) / mean(CVs_max), .groups = "drop")

ref_lines_df <- bind_rows(
  max_ratio_by_dataset %>% mutate(slope = 1,               label = "CVt = CVs")            %>% select(Dataset, slope, label),
  max_ratio_by_dataset %>% mutate(slope = max_ratio_slope, label = "Equal fraction of max") %>% select(Dataset, slope, label)
) %>%
  mutate(intercept = 0)

p_biplot_combined <- ggplot(cv_plot_df, aes(x = mean_CVs, y = mean_CVt,
                                            color = Constituent, shape = Catchment)) +
  geom_abline(data = ref_lines_df, aes(slope = slope, intercept = intercept, linetype = label),
              color = "grey70", linewidth = 0.6, inherit.aes = FALSE) +
  scale_linetype_manual(name = "Reference line",
                        values = c("CVt = CVs" = "dotted", "Equal fraction of max" = "dashed")) +
  geom_errorbar(aes(ymin = mean_CVt - sd_CVt, ymax = mean_CVt + sd_CVt), width = 0) +
  geom_errorbarh(aes(xmin = mean_CVs - sd_CVs, xmax = mean_CVs + sd_CVs), height = 0) +
  geom_point(size = 3) +
  coord_cartesian(xlim = c(0, axis_limit), ylim = c(0, axis_limit)) +
  facet_wrap(~ Dataset) +
  labs(x = "Mean CVs", y = "Mean CVt", color = NULL, shape = NULL,
       title = "CVt vs CVs: real vs. extended data") +
  theme_bw()

p_biplot_combined

ggsave("plots/cvtcvs_extended_means_plot.png", width = 7, height = 3, dpi = 300)



##Time series: real-only vs. extended, combined into one multipanel plot

ts_extended_nm <- real_sites_extended_nm %>%
  group_by(Catchment, Constituent, CampaignNum, Day) %>%
  summarize(Mean_Conc = mean(Conc), Min_Conc = min(Conc), Max_Conc = max(Conc), .groups = "drop") %>%
  mutate(Dataset = "Extended (real + synthetic)")

ts_real_nm <- real_sites_extended_nm %>%
  filter(Source == "real") %>%
  group_by(Catchment, Constituent, CampaignNum, Day) %>%
  summarize(Mean_Conc = mean(Conc), Min_Conc = min(Conc), Max_Conc = max(Conc), .groups = "drop") %>%
  mutate(Dataset = "Real only")

ts_extended_br <- real_sites_extended_br %>%
  group_by(Catchment, Constituent, CampaignNum, Day) %>%
  summarize(Mean_Conc = mean(Conc), Min_Conc = min(Conc), Max_Conc = max(Conc), .groups = "drop") %>%
  mutate(Dataset = "Extended (real + synthetic)")

ts_real_br <- real_sites_extended_br %>%
  filter(Source == "real") %>%
  group_by(Catchment, Constituent, CampaignNum, Day) %>%
  summarize(Mean_Conc = mean(Conc), Min_Conc = min(Conc), Max_Conc = max(Conc), .groups = "drop") %>%
  mutate(Dataset = "Real only")

ts_combined_all <- bind_rows(ts_extended_nm, ts_real_nm, ts_extended_br, ts_real_br)

p_ts_all <- ggplot(ts_combined_all, aes(x = Day / 365.25, y = Mean_Conc, color = Dataset)) +
  geom_ribbon(aes(ymin = Min_Conc, ymax = Max_Conc, fill = Dataset, group = Dataset),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Extended (real + synthetic)" = "grey40",
                                "Real only" = "#2a78d6")) +
  scale_fill_manual(values = c("Extended (real + synthetic)" = "grey40",
                               "Real only" = "#2a78d6")) +
  facet_grid(Constituent ~ Catchment, scales = "free_y") +
  labs(x = "Synthetic year", y = "Concentration across real sites (mean, min-max range)",
       title = "Real-only vs. extended time series", color = NULL, fill = NULL) +
  theme_bw()

p_ts_all
ggsave("plots/cvtcvs_extended_ts_plot.png", width = 10, height = 5, dpi = 300)
