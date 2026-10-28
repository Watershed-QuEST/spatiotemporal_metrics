## Test of variance collapse helper functions
## on toy data sets
# JRB using helper functions from Bre and Sophie

## Import Packages - install any that are FALSE
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse","changepoint"), require, character.only=T)

############################
#### Import & Prep Data ####
############################

# Read in csv - see other files for how to download from google drive
df <- read.csv("drivedata/toy.csv")
df$Date <- as.Date(df$Date, format = "%m/%d/%y")

# Each synoptic "campaign" here spans two consecutive sampling days
# (different subsets of sites visited on each day), so campaigns are
# grouped by sampling month. Adjust this grouping if your data differ.
df$Campaign <- format(df$Date, "%Y-%m")

# convert subcatchment area to km^2 (source data is in m^2)
df$Area.km2 <- df$Area.m2 / 1e6

# create a list of watershed project-campaign specific dataframes
df$Proj_Month <- paste(df$Project,df$Campaign,sep = "-")
l <- split(df, df$Proj_Month)

test <- l$`br-2025-02`

##############################
## Test with helper scripts ##
##############################

## start with breakpoints data prep
source("variance_breakpoints_helper.R")
source("variance_collapse_fx_helper.R")

test_NPOC <- prepare_variance_observations(concentration = test$NPOC..mg.C.L., area = test$Area.km2)
detect_variance_breakpoints(concentration = test$NPOC..mg.C.L., area = test$Area.km2)

#^ might be helpful to also spit out a plot that shows concentration across area as a gut check of whether
# a change point would even be expected
ggplot(test_NPOC, aes(area, scaled_concentration))+geom_point()+theme_bw() # no changepoint detected but kind of seems like there should be?

########################################################################
## Run on data frame with detect_variance_breakpoints_by_group
########################################################################

# First detect variance breakpoints
colnames(df)
NPOC_changepoints <- detect_variance_breakpoints_by_group(data = df, concentration = "NPOC..mg.C.L.", area = "Area.km2", watershed = "Project", event = "Campaign")
TDN_changepoints <- detect_variance_breakpoints_by_group(data = df, concentration = "TDN..mg.N.L.", area = "Area.km2", watershed = "Project", event = "Campaign")

# Then summarize the data
NPOC_collapse <- variance_collapse_by_group(NPOC_changepoints)
TDN_collapse <- variance_collapse_by_group(TDN_changepoints)

#View by event - this is very cool! Nice work
NPOC_collapse$by_event
TDN_collapse$by_event


########################################################################
## PELT / breakpoint sensitivity analysis
########################################################################

# PELT remains the changepoint-search method. Here we vary:
#   1. penalty: AIC, SIC, and MBIC
#   2. minimum observations allowed in each segment: 2 through 5
#
# min_sites is held at 6. Changing min_sites alone only changes whether a
# group is eligible; it does not test sensitivity to removing sites.

sensitivity_settings <- tidyr::crossing(
  penalty = c("AIC", "SIC", "MBIC"),
  min_segment_length = 2:5
)

concentration_columns <- c(
  NPOC = "NPOC..mg.C.L.",
  TDN = "TDN..mg.N.L."
)

run_sensitivity_scenario <- function(penalty, min_segment_length) {
  
  breakpoints <- detect_variance_breakpoints_by_group(
    data = df,
    concentration = concentration_columns,
    area = "Area.km2",
    watershed = "Project",
    event = "Campaign",
    penalty = penalty,
    pen_value = 0,
    min_segment_length = min_segment_length,
    min_sites = 6L
  )
  
  collapse <- variance_collapse_by_group(breakpoints)
  
  collapse$by_event |>
    mutate(
      penalty = penalty,
      min_segment_length = min_segment_length,
      scenario = paste0(
        penalty,
        "; minimum segment length = ",
        min_segment_length
      )
    )
}

# Run every penalty × minimum-segment-length combination
sensitivity_results <- purrr::pmap_dfr(
  sensitivity_settings,
  run_sensitivity_scenario
)

# Distinguish the different reasons a collapse was not reported
sensitivity_results <- sensitivity_results |>
  mutate(
    outcome = case_when(
      detection_status == "too_few_sites" ~
        "Ineligible: too few sites",
      
      detection_status == "zero_concentration_variance" ~
        "Ineligible: zero concentration variance",
      
      collapse_detected ~
        "Variance collapse",
      
      detection_status == "no_changepoint" ~
        "No changepoint",
      
      detection_status == "breakpoints_detected" & !collapse_detected ~
        "Breakpoint without variance decrease",
      
      TRUE ~ status
    ),
    scenario = factor(
      scenario,
      levels = paste0(
        rep(c("AIC", "SIC", "MBIC"), each = 4),
        "; minimum segment length = ",
        rep(2:5, times = 3)
      )
    )
  )

# Inspect the full scenario results
sensitivity_results |>
  select(
    watershed,
    event,
    constituent,
    penalty,
    min_segment_length,
    n_sites,
    detection_status,
    collapse_detected,
    threshold,
    threshold_next_area,
    threshold_midpoint,
    outcome
  ) |>
  arrange(constituent, watershed, event, penalty, min_segment_length)

########################################################################
## Summarize stability for each watershed-event-constituent group
########################################################################

sensitivity_stability <- sensitivity_results |>
  group_by(watershed, event, constituent) |>
  summarize(
    n_scenarios = n(),
    n_collapse = sum(collapse_detected, na.rm = TRUE),
    collapse_fraction = mean(collapse_detected, na.rm = TRUE),
    
    median_threshold = if (any(collapse_detected)) {
      median(threshold[collapse_detected], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    minimum_threshold = if (any(collapse_detected)) {
      min(threshold[collapse_detected], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    maximum_threshold = if (any(collapse_detected)) {
      max(threshold[collapse_detected], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    stability_class = case_when(
      n_collapse == 0 ~ "Never detected",
      n_collapse == n_scenarios ~ "Always detected",
      TRUE ~ "Parameter-sensitive"
    ),
    .groups = "drop"
  )

sensitivity_stability

########################################################################
## Plot 1: overall sensitivity to penalty and segment length
########################################################################

sensitivity_detection_summary <- sensitivity_results |>
  group_by(constituent, penalty, min_segment_length) |>
  summarize(
    collapse_fraction = mean(collapse_detected, na.rm = TRUE),
    .groups = "drop"
  )

p_detection_sensitivity <- ggplot(
  sensitivity_detection_summary,
  aes(
    x = factor(min_segment_length),
    y = penalty,
    fill = collapse_fraction
  )
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("%.2f", collapse_fraction)),
    size = 4
  ) +
  facet_wrap(~constituent) +
  scale_fill_gradient(
    low = "white",
    high = "#2166ac",
    limits = c(0, 1)
  ) +
  labs(
    title = "Sensitivity of variance-collapse detection",
    subtitle = paste(
      "Fraction of watershed-campaign groups classified",
      "as a variance collapse"
    ),
    x = "Minimum observations per segment",
    y = "Penalty",
    fill = "Collapse\nfraction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold")
  )

p_detection_sensitivity

#MBIC is most conservative, and package recommended - continue using that?
#Consistent across min obs per segment, except 5 obs for TDN.

########################################################################
## Plot 2: outcome under every scenario
########################################################################

p_scenario_agreement <- sensitivity_results |>
  mutate(
    group = paste(watershed, event, sep = " / ")
  ) |>
  ggplot(
    aes(
      x = scenario,
      y = group,
      fill = outcome
    )
  ) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(
    ~constituent,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Variance collapse" = "#1b9e77",
      "Breakpoint without variance decrease" = "#d95f02",
      "No changepoint" = "#7570b3",
      "Ineligible: too few sites" = "#bdbdbd",
      "Ineligible: zero concentration variance" = "#636363"
    )
  ) +
  labs(
    title = "Agreement among PELT sensitivity scenarios",
    x = NULL,
    y = "Watershed / campaign",
    fill = "Outcome"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 55,
      hjust = 1
    ),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_scenario_agreement
# I think this is my favorite plot... 


########################################################################
## Example: Detailed diagnostic for BR February 2025 NPOC
########################################################################

example_sensitivity <- sensitivity_results |>
  filter(
    tolower(watershed) == "br",
    event == "2025-02",
    constituent == "NPOC"
  )

# Repeat the prepared observations once for each sensitivity scenario so
# each panel shows the same raw pattern with its scenario-specific threshold.
example_plot_data <- tidyr::crossing(
  test_NPOC,
  sensitivity_settings
)

p_example_sensitivity <- ggplot(
  example_plot_data,
  aes(
    x = area,
    y = scaled_concentration
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey70"
  ) +
  geom_point(
    color = "#2c7fb8",
    size = 2.5
  ) +
  geom_vline(
    data = example_sensitivity |>
      filter(collapse_detected),
    aes(xintercept = threshold),
    color = "#d7301f",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_x_log10() +
  facet_grid(
    penalty ~ min_segment_length,
    labeller = label_both
  ) +
  labs(
    title = "BR February 2025 NPOC: PELT sensitivity",
    subtitle = paste(
      "Dashed red line indicates a detected",
      "variance-collapse threshold"
    ),
    x = expression(Contributing~area~(km^2)),
    y = "Standardized NPOC concentration"
  ) +
  theme_bw(base_size = 11)

p_example_sensitivity
#still leaves which segment length to use unresolved... 

########################################################################
## Optional: save tables and figures
########################################################################

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

write.csv(
  sensitivity_results,
  "outputs/variance_collapse_sensitivity_scenarios.csv",
  row.names = FALSE
)

write.csv(
  sensitivity_stability,
  "outputs/variance_collapse_sensitivity_stability.csv",
  row.names = FALSE
)

ggsave(
  "outputs/variance_collapse_detection_sensitivity.png",
  p_detection_sensitivity,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "outputs/variance_collapse_scenario_agreement.png",
  p_scenario_agreement,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "outputs/variance_collapse_threshold_sensitivity.png",
  p_threshold_sensitivity,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "outputs/BR_2025-02_NPOC_sensitivity.png",
  p_example_sensitivity,
  width = 12,
  height = 9,
  dpi = 300
)
