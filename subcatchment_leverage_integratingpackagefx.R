#Title: Subcatchment leverage#
#This script is for the subcatchment leverage metric to be used for##
##the spatiotemporal metrics collab paper

# Purpose: to plot leverage data from NM and BR for the collab. paper using the semi-toy data
#Author: Andrew Ali, Bre Rivera Waterman
#rerwitten: 08/13/2026 (from prior metric codes) -> package helper fxs implemented 09/24/2026
#last modified: 09/24/2026

#### Packages ####
library(googledrive) 
library(googlesheets4)
library(lubridate)
library(dplyr)
library(tidyverse)
library(ggplot2)

source("subcatchment-leverage-fxhelper.R")


#### Configure in/outputs and file structure ####

data_out_dir <- "data"
plot_dir     <- "plots"

# Constituents to include in this analysis -- must match column names in data
constituents <- c("NPOC..mg.C.L.", "TDN..mg.N.L.")

toy            <- read_csv(file.path(data_out_dir, "NM-BR_Toy_dataset.csv"), show_col_types = FALSE) %>% mutate(Date = mdy(Date))

#preparing data
toy_long <- toy %>%
  pivot_longer(
    cols = all_of(constituents),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    # Use CampaignID here instead if it is the true unique event identifier.
    event = Date,
    month = month(Date, label = TRUE, abbr = TRUE),
    Area.km = Area.m2 / 1e6,
    observation_id = row_number()
  ) %>%
  group_by(Project, event, variable) %>%
  mutate(
    is_outlet = !is.na(Area.m2) &
      Area.m2 == max(Area.m2, na.rm = TRUE)
  ) %>%
  ungroup()

#check maximum-area rule identifies exactly one outlet per group
outlet_check <- toy_long %>%
  summarise(
    n_rows = n(),
    n_outlets = sum(is_outlet),
    outlet_sites = paste(
      unique(as.character(Site[is_outlet])),
      collapse = ", "
    ),
    .by = c(Project, event, variable)
  )

print(outlet_check, n = Inf)

problem_outlets <- outlet_check %>%
  filter(n_outlets != 1)

if (nrow(problem_outlets) > 0) {
  print(problem_outlets, n = Inf)
  
  stop(
    "Every Project-event-variable group must have exactly one maximum-area outlet."
  )
}


#original + dynamic outlet identifier
toy_long <- toy_long %>%
  group_by(Project, event, variable) %>%
  mutate(
    outlet_Area = Area.m2[is_outlet],
    outlet_value = value[is_outlet],
    outlet_Q = Q[is_outlet],
    
    leverage_site = if_else(
      is_outlet,
      ".__dynamic_outlet__.",
      as.character(Site)
    )
  ) %>%
  ungroup()


#### PART A.1 -- Original L + variants calculations ####
leverage_original <- toy_long %>%
  mutate(
    L = (value - outlet_value) *
      (Area.m2 / outlet_Area) *
      (Q / outlet_Q),
    
    Lpct = 100 * L / outlet_value,
    
    LS = (value - outlet_value) *
      (Area.m2 / outlet_Area),
    
    LSpct = 100 * LS / outlet_value
  ) %>%
  select(
    observation_id,
    Project,
    event,
    month,
    variable,
    Site,
    Area.m2,
    Q,
    is_outlet,
    L,
    Lpct,
    LS,
    LSpct
  ) %>%
  pivot_longer(
    cols = c(L, Lpct, LS, LSpct),
    names_to = "leverage_variant",
    values_to = "leverage_original"
  )

print(leverage_original, n = Inf)



#### PART A.2 -- Package fx helper L + variants calculations ####
leverage_package <- map_dfr(
  c("L", "Lpct", "LS", "LSpct"),
  function(current_variant) {
    subcatchment_leverage(
      data = toy_long,
      concentration_col = "value",
      site_col = "leverage_site",
      event_col = "event",
      watershed_col = "Project",
      constituent_col = "variable",
      contributing_area_col = "Area.m2",
      outlet_site = ".__dynamic_outlet__.",
      variant = current_variant,
      specific_discharge_col = if (
        current_variant %in% c("L", "Lpct")
      ) {
        "Q"
      } else {
        NULL
      }
    )
  }
) %>%
  transmute(
    observation_id,
    Project,
    event,
    month,
    variable,
    Site,
    Area.m2,
    Q,
    is_outlet,
    leverage_variant,
    leverage_package = leverage
  )

print(leverage_package, n = Inf)

#### PART A.3 -- Compare original and package results ####
#direct comparison
leverage_comparison <- leverage_original %>%
  full_join(
    leverage_package,
    by = c(
      "observation_id",
      "Project",
      "event",
      "month",
      "variable",
      "Site",
      "Area.m2",
      "Q",
      "is_outlet",
      "leverage_variant"
    )
  ) %>%
  mutate(
    difference = leverage_package - leverage_original,
    absolute_difference = abs(difference),
    agrees = case_when(
      is.na(leverage_original) & is.na(leverage_package) ~ TRUE,
      is.na(leverage_original) | is.na(leverage_package) ~ FALSE,
      TRUE ~ near(
        leverage_original,
        leverage_package,
        tol = 1e-10
      )
    )
  ) %>%
  arrange(
    leverage_variant,
    Project,
    variable,
    event,
    Site
  )

print(leverage_comparison, n = Inf)




#summary x variant
comparison_summary <- leverage_comparison %>%
  summarise(
    n_rows = n(),
    n_compared = sum(
      !is.na(leverage_original) &
        !is.na(leverage_package)
    ),
    n_agree = sum(agrees),
    n_disagree = sum(!agrees),
    maximum_absolute_difference = if (
      all(is.na(absolute_difference))
    ) {
      NA_real_
    } else {
      max(absolute_difference, na.rm = TRUE)
    },
    mean_absolute_difference = if (
      all(is.na(absolute_difference))
    ) {
      NA_real_
    } else {
      mean(absolute_difference, na.rm = TRUE)
    },
    .by = leverage_variant
  )

print(comparison_summary, n = Inf)


#original plots
leverage_plot_data <- leverage_comparison %>%
  mutate(
    Area.km = Area.m2 / 1e6
  ) %>%
  pivot_longer(
    cols = c(leverage_original, leverage_package),
    names_to = "Approach",
    values_to = "leverage"
  ) %>%
  mutate(
    Approach = recode(
      Approach,
      leverage_original = "Original calculation",
      leverage_package = "Package function"
    ),
    Approach = factor(
      Approach,
      levels = c(
        "Original calculation",
        "Package function"
      )
    )
  )


p1 <- leverage_plot_data %>%
  filter(
    leverage_variant == "Lpct",
    !is.na(leverage)
  ) %>%
  ggplot(
    aes(
      x = Area.km,
      y = leverage,
      color = Approach,
      shape = Approach
    )
  ) +
  geom_point(
    alpha = 0.8,
    size = 3.5,
    stroke = 1.1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_grid(
    Project ~ variable,
    scales = "free_y"
  ) +
  scale_x_log10() +
  scale_y_reverse() +
  scale_color_manual(
    values = c(
      "Original calculation" = "#4C78A8",
      "Package function" = "#F58518"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Original calculation" = 1,
      "Package function" = 4
    )
  ) +
  labs(
    title = "Subcatchment leverage by contributing area",
    subtitle = "Comparison of original and package calculations",
    x = "Upstream contributing area (km²)",
    y = "Subcatchment leverage (% of outlet concentration)",
    color = "Calculation approach",
    shape = "Calculation approach"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p1


p2 <- leverage_plot_data %>%
  filter(
    leverage_variant == "LS",
    !is.na(leverage)
  ) %>%
  ggplot(
    aes(
      x = month,
      y = leverage,
      fill = Approach,
      color = Approach
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    alpha = 0.35,
    outlier.shape = NA
  ) +
  geom_point(
    aes(shape = Approach),
    position = position_jitterdodge(
      jitter.width = 0.12,
      jitter.height = 0,
      dodge.width = 0.8,
      seed = 123
    ),
    alpha = 0.65,
    size = 1.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_grid(
    Project ~ variable,
    scales = "free_y"
  ) +
  scale_y_reverse() +
  scale_color_manual(
    values = c(
      "Original calculation" = "#4C78A8",
      "Package function" = "#F58518"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Original calculation" = "#4C78A8",
      "Package function" = "#F58518"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Original calculation" = 16,
      "Package function" = 17
    )
  ) +
  labs(
    title = "Monthly standardized subcatchment leverage",
    subtitle = "Comparison of original and package calculations",
    x = "Month",
    y = "Standardized subcatchment leverage",
    fill = "Calculation approach",
    color = "Calculation approach",
    shape = "Calculation approach"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p2
