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

#############################################################
## Test with helper script - variance_breakpoints_helper.R ##
#############################################################

## start with breakpoints data prep
source("variance_breakpoints_helper.R")

test_NPOC <- prepare_variance_observations(concentration = test$NPOC..mg.C.L., area = test$Area.km2)
detect_variance_breakpoints(concentration = test$NPOC..mg.C.L., area = test$Area.km2)

## Run on data frame with detect_variance_breakpoints_by_group
## grouping by Proj_Month
colnames(df)
NPOC_changepoints <- detect_variance_breakpoints_by_group(data = df, concentration = "NPOC..mg.C.L.", area = "Area.km2", watershed = "Project", event = "Campaign")
TDN_changepoints <- detect_variance_breakpoints_by_group(data = df, concentration = "TDN..mg.N.L.", area = "Area.km2", watershed = "Project", event = "Campaign")

#^ might be helpful to provide guidance that concentration/area/watershed/event should be provided as column names in quotations
#^ also might be helpful to include all information, but then provide a summary table of the key information (i.e., existence of changepoint detected and magnitude if so)









