#This script is used to calculate coefficients of variation with a toy dataset for QuEST's planned commentary manuscript on the use of spatiotemporal metrics for assessing the spatiotemporal variance of stream chemistry across stream networks.
#By: Lauren Giggy (Contact: giggylauren@gmail.com)
#Last Updated: 7/24/26

# CVs --------------------------------------------------------------------
#libraries
library(dplyr)
library(googledrive)
library(lubridate)

### Read in Data
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM_BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")

#Adjust datetime
View(toy)
toy$Date <- lubridate::mdy(toy$Date)

#add area normalized streamflow values - Q must be in units of L/s and area in units of m2 for this formula
toy <- mutate(toy, q_mm_day = (Q* 86400)/Area.m2)

#Subset Toy Data into Projects 
NM_toy <- subset(toy, Project %in% c('nm'))
BR_toy <- subset(toy, Project %in% c('br'))

#Add event IDs to data based on date of data collection. 
#This allows samples collected within 2 days of each other to be combined into 1 event. 
#Change the "units" or "time dif" as necessary for study design. 
NM_toy <-  NM_toy %>%
  arrange(Date) %>%  # make sure data are ordered
  mutate(
    time_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
    new_event = ifelse(is.na(time_diff) | time_diff > 2, 1, 0),
    event_id = cumsum(new_event)
  ) %>%
  dplyr::select(-time_diff, -new_event)

BR_toy <-  BR_toy %>%
  arrange(Date) %>%  # make sure data are ordered
  mutate(
    time_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
    new_event = ifelse(is.na(time_diff) | time_diff > 2, 1, 0),
    event_id = cumsum(new_event)
  ) %>%
  dplyr::select(-time_diff, -new_event)

#Calculate CVs for each synoptic event ID, exclude events with limited sites sampled
#NM
NM_CVs_toy <- NM_toy %>%
  group_by(event_id) %>%
  filter(n() >= 8) %>%   # require at least 8 samples per event_id
  summarise(CVs_TDN = sd(TDN..mg.N.L., na.rm = TRUE) / mean(TDN..mg.N.L., na.rm = TRUE),
            CVs_DOC = sd(NPOC..mg.C.L., na.rm = TRUE) / mean(NPOC..mg.C.L., na.rm = TRUE),
            CVs_q = sd(q_mm_day, na.rm = TRUE) / mean(q_mm_day, na.rm = TRUE),
            survey_start_date = min(Date, na.rm = TRUE), # add dates 
            survey_end_date = max(Date, na.rm = TRUE)) # add dates 
View(NM_CVs_toy)

#BR
BR_CVs_toy <- BR_toy %>%
  group_by(event_id) %>%
  filter(n() >= 8) %>%   # require at least 8 samples per event_id
  summarise(CVs_TDN = sd(TDN..mg.N.L., na.rm = TRUE) / mean(TDN..mg.N.L., na.rm = TRUE),
            CVs_DOC = sd(NPOC..mg.C.L., na.rm = TRUE) / mean(NPOC..mg.C.L., na.rm = TRUE),
            CVs_q = sd(q_mm_day, na.rm = TRUE) / mean(q_mm_day, na.rm = TRUE),
            survey_start_date = min(Date, na.rm = TRUE), # add dates 
            survey_end_date = max(Date, na.rm = TRUE)) # add dates 

View(BR_CVs_toy)



# CVt --------------------------------------------------------------------
#libraries
library(dplyr)
library(googledrive)
library(lubridate)
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM_BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")
View(toy)

#Adjust datetime
View(toy)
toy$Date <- lubridate::mdy(toy$Date)

#add area normalized streamflow values - Q must be in units of L/s and area in units of m2 for this formula
toy <- mutate(toy, q_mm_day = (Q* 86400)/Area.m2)

#Subset Toy Data into Projects 
NM_toy <- subset(toy, Project %in% c('nm'))
BR_toy <- subset(toy, Project %in% c('br'))

#Add event IDs to data based on date of data collection. 
#This allows samples collected within 2 days of each other to be combined into 1 event. 
#Change the "units" or "time dif" as necessary for study design. 
NM_toy <-  NM_toy %>%
  arrange(Date) %>%  # make sure data are ordered
  mutate(
    time_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
    new_event = ifelse(is.na(time_diff) | time_diff > 2, 1, 0),
    event_id = cumsum(new_event)
  ) %>%
  dplyr::select(-time_diff, -new_event)

BR_toy <-  BR_toy %>%
  arrange(Date) %>%  # make sure data are ordered
  mutate(
    time_diff = as.numeric(difftime(Date, lag(Date), units = "days")),
    new_event = ifelse(is.na(time_diff) | time_diff > 2, 1, 0),
    event_id = cumsum(new_event)
  ) %>%
  dplyr::select(-time_diff, -new_event)

#Calculate CVt for each synoptic event ID, exclude sites with limited samples. 
NM_CVt_ToyData <- NM_ToyData %>%
  group_by(Site) %>%
  filter(n() >= 3) %>%   # require at least 3 samples per site
  summarise(CVt_TDN = sd(TDN..mg.N.L., na.rm = TRUE) / mean(TDN..mg.N.L., na.rm = TRUE),
            CVt_DOC = sd(NPOC..mg.C.L., na.rm = TRUE) / mean(NPOC..mg.C.L., na.rm = TRUE),
            CVt_q = sd(q_mm_day, na.rm = TRUE) / mean(q_mm_day, na.rm = TRUE)) 
View(NM_CVt_ToyData)


BR_CVt_ToyData <- BR_ToyData %>%
  group_by(Site) %>%
  filter(n() >= 3) %>%   # require at least 3 samples per site
  summarise(CVt_TDN = sd(TDN..mg.N.L., na.rm = TRUE) / mean(TDN..mg.N.L., na.rm = TRUE),
            CVt_DOC = sd(NPOC..mg.C.L., na.rm = TRUE) / mean(NPOC..mg.C.L., na.rm = TRUE),
            CVt_q = sd(q_mm_day, na.rm = TRUE) / mean(q_mm_day, na.rm = TRUE)) 
View(BR_CVt_ToyData)
