## Variance collapse code demonstration
## Drafted by J.R. Blaszczak with some assistance from Claude

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

########################################
#### Variance collapse calculations ####
########################################

## First, scale the data (z-score) using the whole-watershed mean/SD for this campaign,
# per Shogren et al. (2022): "scaled by subtracting the whole watershed
# mean and dividing by the standard deviation"

f_scale <- function(x){
  
  x$NPOC_z <- as.numeric(scale(x$NPOC..mg.C.L.))
  x$TDN_z <- as.numeric(scale(x$TDN..mg.N.L.))
  return(x)
  
}

l_z <- lapply(l, function(x) f_scale(x))

## Next, order each data frame by ascending subcatchment area
# then use changepoint package to determine changepoint

t <- l_z$`br-2025-03`

t_order <- t[order(t$Area.km2), ]

ggplot(t_order, aes(Area.km2, NPOC_z))+geom_point()


# PELT changepoint in variance, data already centered on 0
# MBIC penalty for PELT (Killick & Eckley default choice)
# minseglen 2 minimum points required on either side of a break
NPOC_fit <- cpt.var(t_order$NPOC_z, method = "PELT", penalty = "MBIC",
               know.mean = TRUE, mu = 0)

cpts(NPOC_fit)

# **having issues with cpts function -- PAUSING HERE - see VarianceCollapse_helper_fxn_test.R




