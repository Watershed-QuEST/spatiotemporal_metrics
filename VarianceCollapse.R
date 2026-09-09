## Variance collapse code demonstration
## Drafted by J.R. Blaszczak

## Import Packages


#### Imports ####
# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM-BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toydata = read.csv("drivedata/toy.csv")

### Wrangle toy data ####
## Group by project
nm_toy <- toydata[which(toydata$Project=='nm'),]

br_toy <- toydata[which(toydata$Project=='br'),]

