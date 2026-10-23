#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-10-07 
# Last update (Person, Date): Alex Webster, 2026-10-08

# This script makes plots for overlaying over Arial's conceptual figure for integrated interpretation of multiple metrics. 
# See summary of values from toy datasets in this Google Sheet: https://docs.google.com/spreadsheets/d/1GgaCtrxUyJ9JjP1zAec6U1-qWqFGc6Z0CwxlNxzvbXQ/edit?usp=sharing

# Outputs:
#   plots/*_means_plot.png with * denoting the metric comparison

#### Packages ####
library(tidyverse)

#### CVs and CVt ####

# summary stats using just the real toy datasets from br and nm, where CVt samples size = 3 and CVs sample size = ~20
cv_summary <- tribble(
  ~Metric,    ~Constituent, ~CVs,   ~CVt,
  "Mean BR",  "DOC",        0.219,  0.144,
  "SD BR",    "DOC",        0.042,  0.0602,
  "Mean NM",  "DOC",        0.501,  0.209,
  "SD NM",    "DOC",        0.114,  0.229,
  "Mean BR",  "TDN",        0.517,  0.211,
  "SD BR",    "TDN",        0.035,  0.156,
  "Mean NM",  "TDN",        0.724,  0.244,
  "SD NM",    "TDN",        0.114,  0.124
)


# summary stats using the real toy datasets + synthetic data at matching sites to extend the dataset to 20 consecutive campaigns. We are doing this so that CVt has the same sample size as CVs so that the sample size does not impact the interpretation of whether spatial or temporal variance is stronger. See CVt.R for code that pulled the extended dataset and summarized it. 
cv_extended_summary <- tribble(
  ~Metric,    ~Constituent, ~CVs,   ~CVt,
  "Mean BR",  "DOC",        0.289,  1.22,
  "SD BR",    "DOC",        0.0417, 0.0947,
  "Mean NM",  "DOC",        0.466,  1.13,
  "SD NM",    "DOC",        0.0830, 0.0971,
  "Mean BR",  "TDN",        0.832,  1.01,
  "SD BR",    "TDN",        0.240,  0.179,
  "Mean NM",  "TDN",        0.987,  0.509,
  "SD NM",    "TDN",        0.264,  0.114
)
# note that the CVs (here calculated including the synthetic data) is very similar to that of the real toy dataset. In contrast, CVt is much higher now.

cv_plot_df = cv_summary %>% 
  separate(Metric, into = c("Stat", "Catchment"), sep = " ") %>%
  pivot_longer(c(CVs, CVt), names_to = "Metric", values_to = "value") %>%
  pivot_wider(names_from = c(Stat, Metric), values_from = value)

cv_plot_df_extended <- cv_extended_summary %>% # 
  separate(Metric, into = c("Stat", "Catchment"), sep = " ") %>%
  pivot_longer(c(CVs, CVt), names_to = "Metric", values_to = "value") %>%
  pivot_wider(names_from = c(Stat, Metric), values_from = value)


# calculate axes limits
## ---- Maximum possible CV (sqrt(n), sample SD) ----
cv_max <- function(n) sqrt(n)

## Sample sizes behind each CV, per catchment
##   n_sites     = sites per campaign (sets the CVs ceiling)
##   n_campaigns = campaigns per site (sets the CVt ceiling)
## Replace the site counts with your actual numbers
cv_n <- tribble(
  ~Catchment, ~n_sites, ~n_campaigns,
  "NM",       22,       3,
  "BR",       19,       3
)

cv_n_extended <- tribble(
  ~Catchment, ~n_sites, ~n_campaigns,
  "NM",       22,       20,
  "BR",       19,       20
)

cv_plot_df <- cv_plot_df %>%
  left_join(cv_n, by = "Catchment") %>%
  mutate(CVs_max = cv_max(n_sites),
         CVt_max = cv_max(n_campaigns),
         Mean_CVs_frac = Mean_CVs / CVs_max,
         Mean_CVt_frac = Mean_CVt / CVt_max,
         SD_CVs_frac   = SD_CVs / CVs_max,
         SD_CVt_frac   = SD_CVt / CVt_max)

cv_plot_df

cv_plot_df_extended <- cv_plot_df_extended %>%
  left_join(cv_n_extended, by = "Catchment") %>%
  mutate(CVs_max = cv_max(n_sites),
         CVt_max = cv_max(n_campaigns),
         Mean_CVs_frac = Mean_CVs / CVs_max,
         Mean_CVt_frac = Mean_CVt / CVt_max,
         SD_CVs_frac   = SD_CVs / CVs_max,
         SD_CVt_frac   = SD_CVt / CVt_max)

cv_plot_df_extended


## regular plots for immediate viewing

# real toy data only
ggplot(cv_plot_df, aes(x = Mean_CVs, y = Mean_CVt,
                       shape = Catchment, linetype = Catchment, color=Constituent)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey70", linewidth = 0.6) +   # raw equality
  geom_abline(slope = mean(cv_plot_df$CVt_max) / mean(cv_plot_df$CVs_max), intercept = 0,
              linetype = "dashed", color = "grey70", linewidth = 0.6) +   # equal fraction of max
  geom_errorbar(aes(ymin = Mean_CVt - SD_CVt, ymax = Mean_CVt + SD_CVt),
                width = 0.025, linewidth = 1.2) +
  geom_errorbar(aes(xmin = Mean_CVs - SD_CVs, xmax = Mean_CVs + SD_CVs),
                width = 0.025, linewidth = 1.2,
                orientation = "y") +
  geom_point(size = 4, stroke = 1.5) +
  scale_shape_manual(values = c(BR = 16, NM = 2)) +          # filled circle vs. open triangle
  scale_linetype_manual(values = c(BR = "solid", NM = "22")) + # solid vs. dashed
  coord_equal(xlim = c(-0.05, 1.5), ylim = c(-0.05, 1.5), expand = FALSE) +
  labs(x = "Spatial CV (CVs)",
       y = "Temporal CV (CVt)",
       shape = "Catchment", linetype = "Catchment") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.key.width = unit(1.5, "cm"),
    panel.grid = element_blank()
  )

# real + extended data 
ggplot(cv_plot_df_extended, aes(x = Mean_CVs, y = Mean_CVt,
                       shape = Catchment, linetype = Catchment, color=Constituent)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey70", linewidth = 0.6) +   # raw equality
  geom_abline(slope = mean(cv_plot_df_extended$CVt_max) / mean(cv_plot_df_extended$CVs_max), intercept = 0,
              linetype = "dashed", color = "grey70", linewidth = 0.6) +   # equal fraction of max
  geom_errorbar(aes(ymin = Mean_CVt - SD_CVt, ymax = Mean_CVt + SD_CVt),
                width = 0.025, linewidth = 1.2) +
  geom_errorbar(aes(xmin = Mean_CVs - SD_CVs, xmax = Mean_CVs + SD_CVs),
                width = 0.025, linewidth = 1.2,
                orientation = "y") +
  geom_point(size = 3, stroke = 1.5) +
  scale_shape_manual(values = c(BR = 16, NM = 15)) +          # filled circle vs. open triangle
  #scale_linetype_manual(values = c(BR = "solid", NM = "22")) + # solid vs. dashed
  coord_equal(xlim = c(-0.05, 1.5), ylim = c(-0.05, 1.5), expand = FALSE) +
  labs(x = "Spatial CV (CVs)",
       y = "Temporal CV (CVt)",
       shape = "Catchment", linetype = "Catchment") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.key.width = unit(1.5, "cm"),
    panel.grid = element_blank()
  )



## transparent plot for overlay

# real toy data only
ggplot(cv_plot_df, aes(x = Mean_CVs, y = Mean_CVt,
                       shape = Catchment, linetype = Catchment, color=Constituent)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey70", linewidth = 0.6) +   # raw equality
  geom_abline(slope = mean(cv_plot_df$CVt_max) / mean(cv_plot_df$CVs_max), intercept = 0,
              linetype = "dashed", color = "grey70", linewidth = 0.6) +   # equal fraction of max
  geom_errorbar(aes(ymin = Mean_CVt - SD_CVt, ymax = Mean_CVt + SD_CVt),
                width = 0.025, linewidth = 1.2, color = "white") +
  geom_errorbar(aes(xmin = Mean_CVs - SD_CVs, xmax = Mean_CVs + SD_CVs),
                width = 0.025, linewidth = 1.2, color = "white",
                orientation = "y") +
  geom_point(size = 4, color = "white", fill = "white", stroke = 1.5) +
  scale_shape_manual(values = c(BR = 16, NM = 2)) +          # filled circle vs. open triangle
  scale_linetype_manual(values = c(BR = "solid", NM = "22")) + # solid vs. dashed
  coord_equal(xlim = c(-0.05, 1.5), ylim = c(-0.05, 1.5), expand = FALSE) +
  labs(x = "",
       y = "",
       shape = "Catchment", linetype = "Catchment") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.key.width = unit(1.5, "cm"),
    # white text and lines so they show on a dark background
    text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    axis.ticks = element_line(color = "white"),
    panel.border = element_rect(color = "white", fill = NA),
    panel.grid = element_blank()
  )
ggsave("plots/cvtcvs_real_means_plot.png", width = 3, height = 2.5, dpi = 300, bg = "transparent")


# real + extended data 
ggplot(cv_plot_df_extended, aes(x = Mean_CVs, y = Mean_CVt,
                       shape = Catchment, linetype = Catchment, color=Constituent)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey70", linewidth = 0.6) +   # raw equality
  #geom_abline(slope = mean(cv_plot_df_extended$CVt_max) / mean(cv_plot_df_extended$CVs_max), intercept = 0,
   #           linetype = "dashed", color = "grey70", linewidth = 0.6) +   # equal fraction of max
  geom_errorbar(aes(ymin = Mean_CVt - SD_CVt, ymax = Mean_CVt + SD_CVt),
                width = 0.025, linewidth = 1.2, color = "white") +
  geom_errorbar(aes(xmin = Mean_CVs - SD_CVs, xmax = Mean_CVs + SD_CVs),
                width = 0.025, linewidth = 1.2, color = "white",
                orientation = "y") +
  geom_point(size = 2, color = "white", fill = "white", stroke = 1.5) +
  scale_shape_manual(values = c(BR = 16, NM = 15)) +          # filled circle vs. open triangle
  #scale_linetype_manual(values = c(BR = "solid", NM = "22")) + # solid vs. dashed
  coord_equal(xlim = c(-0.05, 1.5), ylim = c(-0.05, 1.5), expand = FALSE) +
  labs(x = "",
       y = "",
       shape = "Catchment", linetype = "Catchment") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.key.width = unit(1.5, "cm"),
    # white text and lines so they show on a dark background
    text = element_text(color = "white"),
    axis.text = element_text(color = "white"),
    axis.ticks = element_line(color = "white"),
    panel.border = element_rect(color = "white", fill = NA),
    panel.grid = element_blank()
  )
ggsave("plots/cvtcvs_extended_means_plot.png", width = 3, height = 2.5, dpi = 300, bg = "transparent")
