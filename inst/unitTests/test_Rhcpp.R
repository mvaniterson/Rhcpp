library(Rhcpp)

##load data
data(hcpp)

##set parameters
F <- hcpp$F
Y <- hcpp$Y
k <- 10
lambda1 <- 20
lambda2 <- 1
lambda3 <- 1
iter <- 100

##read expected result
Mres <- read.table(file.path(path.package("Rhcpp"), "extdata", "Res.txt.bz2"), sep="\t")
res <- hcp(F, Y, k, lambda1, lambda2, lambda3, iter, fast=FALSE)

library(RUnit)
R <- as.vector(res$Res)
M <- as.numeric(unlist(t(Mres)))
checkEquals(R, M, tolerance=0.001) ##not exact probably due to errors differences
checkEquals(R, M) ##not exact probably due to errors differences

res <- hcp(F, Y, k, lambda1, lambda2, lambda3, iter)
R <- as.vector(res$Res)
checkEquals(R, M, tolerance=0.05) ##slightly higher because of different initialization


