#Title: Subcatchment leverage#
#This script is for the subcatchment leverage metric to be used for##
##the spatiotemporal metrics collab paper

# Purpose: to plot leverage data from NM and BR for the collab. paper using the semi-toy data
#Author: Andrew Ali 
#rerwitten: 08/13/2026 (from prior metric codes)
#last modified: 08/13/2026

##############
## Required ackages ##
library(googledrive) 
library(googlesheets4)
library(lubridate)
library(dplyr)
library(tidyverse)
library(ggplot2)

### Read in Data
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM-BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")

#Adjust datetime
View(toy)
toy$Date <- lubridate::mdy(toy$Date)


#Creating a longer dataframe, columns for outlet_Area.m2 and outlet_value to enable me calculate leverage
toy_long <- toy %>% 
  pivot_longer(
    cols = c(NPOC..mg.C.L., TDN..mg.N.L.),
    names_to = "variable",
    values_to = "value"
  ) %>% 
  mutate(month = month(as.Date(Date), label = TRUE, abbr = TRUE)) %>% 
  group_by(Project, month, variable) %>% 
  mutate(
    outlet_Area = max(Area.m2, na.rm = T),
    outlet_value = value[which.max(Area.m2)],
    outlet_Q = Q[which.max(Area.m2)],
    Area.km = Area.m2/1e6,
    variable = factor(variable, levels = c("NPOC..mg.C.L.", "TDN..mg.N.L.")),
    Project = factor(Project, levels = c("nm", "br"))
  ) %>% 
  ungroup()

##################
#Calculating the different leverage metrics based on Table 1
#legend: L = leverage in conc.; LS = leverage when Spec Q is 1; etc
toy_plot <- toy_long %>% 
  mutate(
    L = ((value - outlet_value) * (Area.m2/outlet_Area) * (Q/outlet_Q)),
    L_pct = (L/outlet_value) * 100,
    LS = ((value - outlet_value) * (Area.m2/outlet_Area) * 1),
    LS_pct = (LS/outlet_value) * 100
  )

#partitioning samples by Project
Lev_nm <- subset(toy_plot, Project == "nm")
Lev_br <- subset(toy_plot, Project == "br")


###############
#plotting leverage 
p1 <- ggplot(toy_plot, aes(x = Area.km, y = L_pct)) + #we can change y axis to L, LS, LS_pct
  #geom_point(aes(shape = variable), alpha = 0.75, size = 5) +
  geom_jitter(aes(shape = variable), width = 0.15, height = 0.02, alpha = 0.75,
    size = 4) +
  geom_hline(yintercept=0, linetype = 'dashed') +
  facet_grid(Project ~ variable, scales = "free_y") +
  scale_x_log10() + 
  labs(
    title = "Subcatchment Leverage by Area (DOC & TDN)",
    x = "Subcatchment area (km²)",
    y = expression("Subcatchment leverage(%)")
  ) +
  scale_y_reverse()+
  #ylim(c(2.5, -1))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )
p1

# Boxplot of leverage by month
#just as with codes above, we can switch the y-axis for the variant of leverage
p2 <- ggplot(toy_plot, aes(x = month, y = LS)) +
  geom_boxplot(fill = "lightblue", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(Project ~ variable, scales = "free_y") +
  scale_y_reverse() +
  labs(
    title = "Monthly net behavior in subcatchment leverage",
    x = "Month",
    y = "Subcatchment leverage"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )
p2
