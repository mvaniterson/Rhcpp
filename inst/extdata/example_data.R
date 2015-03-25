library(Rhcpp)
library(R.matlab)
RhcppDevelPath <- "/media/Storage/packages/Rhcpp"
data <- readMat(file.path(RhcppDevelPath, "inst/extdata", "example_data.mat"))
Y <- data$Y
F <- data$F
save(data, file=file.path(RhcppDevelPath, "data", "hcpp.RData"))
