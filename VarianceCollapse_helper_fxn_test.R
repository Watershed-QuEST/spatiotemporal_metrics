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









