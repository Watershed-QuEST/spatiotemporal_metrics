### load data ###
data_out_dir <- "data"
plot_dir     <- "plots"

n_iter    <- 10000                                  # Monte Carlo iterations
site_grid <- c(3:50, 75, 100, 150, 300)              # site counts to test

actual_n_sites <- 23  # real number of nm sites sampled -- update if that changes

set.seed(42)

clean            <- read_csv(file.path(data_out_dir, "nm_clean.csv"), show_col_types = FALSE)
field_setups     <- readRDS(file.path(data_out_dir, "nm_field_setups.rds"))
synthetic_extended <- read_csv(file.path(data_out_dir, "nm_synthetic_extended.csv"), show_col_types = FALSE)

# ### Amy Leverage equation ### -------------------------------------------

#Determine scaling constant between Q and Area
# 1) create a dataframe that is the mean Q at each site 
mean_flow_NM <- clean %>%
  group_by(Site) %>%
  summarise(mean_Q = mean(Q, na.rm = TRUE), 
            Area.m2 = max(Area.m2, na.rm = TRUE), .groups = 'drop') %>%
  as.data.frame()

# 2) plot in log-log space
ggplot(data = mean_flow_NM, aes(x = log(Area.m2), y = log(mean_Q) ))+
  geom_hline(yintercept = 0)+
  # scale_y_log10()+
  # scale_x_log10()+
  stat_poly_line(se = FALSE) +
  stat_poly_eq(
    aes(label = ifelse(
      after_stat(p.value),
      paste(after_stat(eq.label), after_stat(p.value.label), after_stat(n.label), sep = "*\", \"*"),
      NA)),
    color = "black",  size = 3) +
  geom_point(size =3, shape = 16, alpha = 0.7)+
  theme_classic(base_size =15 )+
  theme(legend.position="right")

#3) Pull b value from log-log relationship
flow_area_model <- lm(log(mean_Q) ~ log(Area.m2), data = mean_flow_NM)

b <- coef(flow_area_model)[["log(Area.m2)"]]
b


# Add nessecary columns to original toy data frame
toyNMdata_leverage <- clean #rename df for clarity

#add Outlet Q column 
toyNMdata_leverage <- toyNMdata_leverage %>%
  left_join(toyNMdata_leverage %>%
      filter(Site == "USF12") %>%
      select(CampaignID, Outlet_Q = Q), by = "CampaignID")

#add Outlet DOC and TDN columns
toyNMdata_leverage <- toyNMdata_leverage %>%
  left_join(toyNMdata_leverage %>%
              filter(Site == "USF12") %>%
              select(CampaignID, Outlet_DOC = NPOC..mg.C.L.), by = "CampaignID")

toyNMdata_leverage <- toyNMdata_leverage %>%
  left_join(toyNMdata_leverage %>%
              filter(Site == "USF12") %>%
              select(CampaignID, Outlet_TDN = TDN..mg.N.L.), by = "CampaignID")
#create variable  = outlet catchment area
NM_Outlet_area <- toyNMdata_leverage %>%
  filter(Site == "USF12") %>%
  pull(Area.m2) %>%
  unique()


# apply Amy's equation to toy dataset
# where "NPOC..mg.C.L." is site Concentration, "Outlet_DOC" is outlet concentration,
# "Area.m2" is site's catchment area, "Outlet_Area.m2" is Outlet's catchment area, 
# and b = scaling constant determined above

toyNMdata_leverage <- mutate(toyNMdata_leverage, DOC_leverage_amy = 
                               (NPOC..mg.C.L.*(Q/ Outlet_Q)) - ((Outlet_DOC)*((Area.m2/NM_Outlet_area)^b)))

toyNMdata_leverage <- mutate(toyNMdata_leverage, TDN_leverage_amy = 
                               (TDN..mg.N.L.*(Q/Outlet_Q)) - ((Outlet_TDN *(Area.m2/NM_Outlet_area)^b)))

#plot results
amy_leverage<- ggplot(data = toyNMdata_leverage, aes(x = Area.m2, y = DOC_leverage_amy, color= CampaignID))+
  geom_hline(yintercept = 0)+
  geom_point(size =3, shape = 16, alpha = 0.7)+
  theme_classic(base_size =15 )+
  theme(legend.position="right")
amy_leverage


# ### Typical Leverage equation ### -----------------------------------------------

#add column that is area-normalized discharge (aka specific discharge, runoff)
#converts L/s to mm/s and then to mm/day for easier sanity check
toyNMdata_leverage <- mutate(toyNMdata_leverage, q_mm_day = ((Q)/Area.m2)*86400) #conversion assumes Q is in L/sec

#add other nessecary columns if not already present from above

#add Outlet q_mm_day column
toyNMdata_leverage <- toyNMdata_leverage %>%
  left_join(toyNMdata_leverage %>%
              filter(Site == "USF12") %>%
              select(CampaignID, Outlet_q_mm_day = q_mm_day), by = "CampaignID")
# 
# #add Outlet DOC and TDN columns
# toyNMdata_leverage <- toyNMdata_leverage %>%
#   left_join(toyNMdata_leverage %>%
#               filter(Site == "USF12") %>%
#               select(CampaignID, Outlet_DOC = NPOC..mg.C.L.), by = "CampaignID")
# 
# toyNMdata_leverage <- toyNMdata_leverage %>%
#   left_join(toyNMdata_leverage %>%
#               filter(Site == "USF12") %>%
#               select(CampaignID, Outlet_TDN = TDN..mg.N.L.), by = "CampaignID")
# #create variable  = outlet catchment area
# NM_Outlet_area <- toyNMdata_leverage %>%
#   filter(Site == "USF12") %>%
#   pull(Area.m2) %>%
#   unique()

#apply typical leverage equation
toyNMdata_leverage <- mutate(toyNMdata_leverage, DOC_leverage = 
                               ((NPOC..mg.C.L.-Outlet_DOC)*(Area.m2/NM_Outlet_area)*(q_mm_day/Outlet_q_mm_day)))
toyNMdata_leverage <- mutate(toyNMdata_leverage, TDN_leverage = 
                               (TDN..mg.N.L.-Outlet_TDN)*(Area.m2/NM_Outlet_area)*(q_mm_day/Outlet_q_mm_day))



# plot results
typical_leverage <- ggplot(data = toyNMdata_leverage, aes(x = Area.m2, y = DOC_leverage, color= CampaignID))+
   geom_hline(yintercept = 0)+
   geom_point(size =3, shape = 16, alpha = 0.7)+
   theme_classic(base_size =15 )+
   theme(legend.position="right")
typical_leverage
##
#compare with Amy's values
plot_grid(typical_leverage, amy_leverage)
