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

# I still think some visualization output options would be helpful - maybe the changepoint detection settings are too sensitive?
