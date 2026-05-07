#This script is used to calculate source area residuals with a toy dataset for QuEST's planned commentary manuscript on the use of spatiotemporal metrics for assessing the spatiotemporal variance of stream chemistry across stream networks.

#### Load Packages ####
library(tidyverse)
library(googledrive)
library(lubridate)
library(dataRetrieval)
library(ggpmisc)

#### Read in data ####

# List all files in the folder
toy_files <- drive_ls(drive_get("https://drive.google.com/drive/u/1/folders/1zh0YTDM5w971iFwmw-iSyTDQQ4MyGL8-"))
# Download the CSV file
googledrive::drive_download(file = toy_files$id[toy_files$name=="NM_BR Toy dataset.csv"], 
                            path = "drivedata/toy.csv",
                            overwrite = T)
# read in csv
toy = read.csv("drivedata/toy.csv")

# add instananeous mass loads
toy$tdn_mass = toy$Q*toy$TDN..mg.N.L.
toy$npoc_mass = toy$Q*toy$NPOC..mg.C.L.

# add watershed ids
toy$ws = NA
toy$ws[1:78] = "NM"
toy$ws[79:137] = "BR"

#### Create log-linear master regressions ####

# --- List of solute mass columns ---
mass_cols <- c("npoc_mass",  "tdn_mass")

# --- Compute site-level mean per solute and Q (for regression) ---
site_summary <- toy %>%
  group_by(Site) %>%
  summarise(
    ws        = first(ws),
    area_m2        = first(Area.m2),
    
    across(all_of(c("Q", mass_cols)),
           ~mean(.x, na.rm = TRUE),
           .names = "{.col}_mean"),
    .groups = "drop"
  ) %>%
  mutate(
    logArea = log10(area_m2),
    
    logQ    = log10(ifelse(Q_mean > 0, Q_mean, NA)),
    logNPOC = log10(ifelse(npoc_mass_mean > 0, npoc_mass_mean, NA)),
    logTDN  = log10(ifelse(tdn_mass_mean > 0, tdn_mass_mean, NA))
  )

# --- Fit log-log linear models ---
lm_list <- list(
  logQ    = lm(logQ ~ logArea, data = site_summary),
  logNPOC = lm(logNPOC ~ logArea, data = site_summary),
  logTDN  = lm(logTDN ~ logArea, data = site_summary)
)

# --- Graph example (DOC) ---
ggplot(site_summary[site_summary$ws=="NM",], aes(x = logArea, y = logNPOC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5
  ) +
  labs(
    x = expression(Log[10]~Subcatchment~Area~(m^2)),
    y = expression(Log[10]~Site~Avg.~DOC~Mass),
    title = "NM: Log–Log Relationship Between Site-Averaged DOC Mass \nand Subcatchment Area"
  ) +
  theme_classic(base_size = 14)
ggplot(site_summary[site_summary$ws=="BR",], aes(x = logArea, y = logNPOC)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5
  ) +
  labs(
    x = expression(Log[10]~Subcatchment~Area~(m^2)),
    y = expression(Log[10]~Site~Avg.~DOC~Mass),
    title = "BR: Log–Log Relationship Between Site-Averaged DOC Mass \nand Subcatchment Area"
  ) +
  theme_classic(base_size = 14)

#### calculate residuals for individual synoptic campaigns ####

# --- Add log-transformed values to full dataset ---
toy_log <- toy %>%
  mutate(
    logArea = log10(Area.m2),
    logQ    = log10(ifelse(Q > 0, Q, NA)),
    logNPOC = log10(ifelse(npoc_mass > 0, npoc_mass, NA)),
    logTDN  = log10(ifelse(tdn_mass > 0, tdn_mass, NA)
  ))

# --- Predict values + residuals ---
for(solute in names(lm_list)) {
  
  pred_name <- paste0("PL_", solute)
  res_name  <- paste0("res_", solute)
  
  toy_log[[pred_name]] <- predict(lm_list[[solute]], newdata = toy_log)
  
  toy_log[[res_name]] <-
    toy_log[[solute]] - toy_log[[pred_name]]
}

# --- Summarize residuals per site across campaigns ---
res_cols <- grep("^res_", names(toy_log), value = TRUE)

residual_summary <- toy_log %>%
  group_by(Site, ws, Area.m2) %>%
  summarise(
    across(all_of(res_cols),
           list(
             mean = ~mean(.x, na.rm = TRUE),
             se   = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
           ),
           .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# --- Plot residuals ---

ggplot(residual_summary[residual_summary$ws=="NM",])+
  geom_bar( aes(x=reorder(Site, res_logTDN_mean), y=res_logTDN_mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Site, ymin=res_logTDN_mean-res_logTDN_se, ymax=res_logTDN_mean+res_logTDN_se), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  ylim(-3, 1)

ggplot(residual_summary[residual_summary$ws=="BR",])+
  geom_bar( aes(x=reorder(Site, res_logTDN_mean), y=res_logTDN_mean), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=Site, ymin=res_logTDN_mean-res_logTDN_se, ymax=res_logTDN_mean+res_logTDN_se), width=0.4, colour="orange", alpha=0.9, linewidth=1.3)+
  ylim(-3, 1.2)
