#### READ ME ####

# Project: Spatiotemporal Metrics
# Author: Alex Webster, 2026-10-07 
# Last update (Person, Date): Alex Webster, 2026-10-08

# This script makes plots for overlaying over Ariel's conceptual figure for integrated interpretation of multiple metrics. 
# See summary of values from toy datasets in this Google Sheet: https://docs.google.com/spreadsheets/d/1GgaCtrxUyJ9JjP1zAec6U1-qWqFGc6Z0CwxlNxzvbXQ/edit?usp=sharing

# Outputs:
#   plots/*_means_plot.png with * denoting the metric comparison

#### Packages ####
library(tidyverse)

#### CVs and CVt ####

cv_summary <- tribble(
  ~Metric,        ~CVs,   ~CVt,
  "Mean DOC NM",  0.501,  0.209,
  "Mean DOC BR",  0.219,  0.144,
  "SD DOC BR",    0.0424, 0.0602,
  "SD DOC NM",    0.114,  0.229
)

cv_extended_summary <- tribble(
  ~Metric,        ~CVs,   ~CVt,
  "Mean DOC NM",  0.466,  1.13,
  "Mean DOC BR",  0.289,  1.22,
  "SD DOC BR",    0.0417, 0.0947,
  "SD DOC NM",    0.0830,  0.0971
)

cv_plot_df <- cv_extended_summary %>% # or put in cv_extended_summary
  separate(Metric, into = c("Stat", "Constituent", "Catchment"), sep = " ") %>%
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

ggplot(cv_plot_df, aes(x = Mean_CVs, y = Mean_CVt,
                       shape = Catchment, linetype = Catchment)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey70", linewidth = 0.6) +   # raw equality
  geom_abline(slope = 1.73 / 4.525, intercept = 0,
              linetype = "dashed", color = "grey70", linewidth = 0.6) +   # equal fraction of max
  geom_errorbar(aes(ymin = Mean_CVt - SD_CVt, ymax = Mean_CVt + SD_CVt),
                width = 0.025, linewidth = 1.2, color = "white") +
  geom_errorbar(aes(xmin = Mean_CVs - SD_CVs, xmax = Mean_CVs + SD_CVs),
                width = 0.025, linewidth = 1.2, color = "white",
                orientation = "y") +
  geom_point(size = 4, color = "white", fill = "white", stroke = 1.5) +
  scale_shape_manual(values = c(BR = 16, NM = 2)) +          # filled circle vs. open triangle
  scale_linetype_manual(values = c(BR = "solid", NM = "22")) + # solid vs. dashed
  coord_equal(xlim = c(-0.05, 2), ylim = c(-0.05, 2), expand = FALSE) +
  labs(x = "Spatial CV (CVs)",
       y = "Temporal CV (CVt)",
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

ggsave("plots/cvtcvs_means_plot.png", width = 3, height = 2.5, dpi = 300, bg = "transparent")
