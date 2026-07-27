#libraries
library(dplyr)
library(googledrive)
library(lubridate)
library(Hmisc)
library(broom)
library(purrr)
library(tidyverse)


### Read in Data
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM_BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")


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
#View(toy)
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

#Check for repeats
NM_toy |>
  summarise(n = dplyr::n(), .by = c(Site, event_id)) |>
  filter(n > 1)

BR_toy |>
  summarise(n = dplyr::n(), .by = c(Site, event_id)) |>
  filter(n > 1)


### median spearman rank and summary table function 
calc_event_cor <- function(data, analyte, min_samples = 8) { 
  #adjust min_samples as preferred
  
  # Create Solute × Event matrix and set up for correlation calculation
  mat <- data %>%
    dplyr::select(Site, event_id, all_of(analyte)) %>% #keep only the columns needed for calculation
    drop_na() %>% #remove NAs
    pivot_wider(names_from = event_id, values_from = all_of(analyte)) %>% #rotate the dataframe so each column is a synoptic event
    arrange(Site) %>% #arrange matching Sites
    select(-Site) %>% #removes Site for correlation 
    as.matrix()
  
  # Calculate the number of paired observations for each event
  n_mat <- outer(
    seq_len(ncol(mat)),
    seq_len(ncol(mat)),
    Vectorize(function(i, j) {
      sum(complete.cases(mat[, c(i, j)]))
    }))
  
  # Calculate event-by-event Spearman correlations
  cor_mat <- cor(mat, method = "spearman",
                 use = "pairwise.complete.obs")
  
  # Remove correlations based on too few paired observations
  cor_mat[n_mat < min_samples] <- NA
  
  #Remove event correlations with itself 
  diag(cor_mat) <- NA
  
  #create a table of the median spearman rank, MARGIN = 2 means it will apply across columns (= 1 applies across rows)
  tibble(event_id = as.numeric(colnames(cor_mat)),
         median_spearman = apply(cor_mat, MARGIN = 2, median, na.rm = TRUE), analyte = analyte)} 

#define analytes
analytes <- c(
  "NPOC..mg.C.L.",
  "TDN..mg.N.L.",
  "q_mm_day")

#Run function on data
NM_toy_sp_pairs <- map_dfr(analytes, ~calc_event_cor(NM_toy, .x))
BR_toy_sp_pairs <- map_dfr(analytes, ~calc_event_cor(BR_toy, .x))

