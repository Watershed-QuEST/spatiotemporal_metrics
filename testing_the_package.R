

#Here is how you can get the package working on your device!!

install.packages("devtools")

library(devtools)

install_github("sjterian/Watershed-metrics") #creator/repo name

library(watershedmetrics)


#Test a function:

spatial_cv <- spatial_cv(c(0.8, 1.1, 1.4, 1.0))
print(spatial_cv)


