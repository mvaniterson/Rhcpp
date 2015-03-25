library(Rhcpp)

##load data
data(hcpp)

##set parameters
F <- data$F
Y <- data$Y
k <- 10
lambda1 <- 20
lambda2 <- 1
lambda3 <- 1
iter <- 100

##read expected result
Mres <- read.table(file.path(path.package("Rhcpp"), "extdata", "Res.txt.bz2"), sep="\t")

Rres <- r_hcp(Y, F, k, lambda1, lambda2, lambda3, iter)
str(Rres)

rrss <- function(obs, exp)
  {
    obs <-  obs$Yn - obs$Z%*%obs$B
    sum((as.vector(Res) - as.vector(exp))^2)
  }

checkEquals(rrss(Rres, Mres), 0)

Rcpparmares <- Rcpparma_hcp(Y, F, k, lambda1, lambda2, lambda3, iter)

checkEquals(rrss(Rcpparmares, Mres), 0)
