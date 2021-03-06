% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hcpcv.R
\name{hcpcv}
\alias{hcpcv}
\title{Perform the HCP normalization algorithm on a grid of model parameters}
\usage{
hcpcv(Z, Y, kRange = c(10, 20), lambdaRange = c(1, 5, 10, 20),
  performance = NULL, iter = 100, stand = TRUE, log = TRUE,
  verbose = TRUE, fast = TRUE)
}
\arguments{
\item{Z}{a matrix nxd of known covariates, where n is the number of
subjects and d is the number of known covariates. *must be standardize
(columns have 0 mean and constant SS).}

\item{Y}{a matrix of nxg of expression data (must be standardized (columns
scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.}

\item{kRange}{multiple numbers of inferred hidden components (k is an integer)}

\item{lambdaRange}{multiple model parameters}

\item{performance}{function accepting res with res$Res the transformed Residuals}

\item{iter}{(optional) iter: number of iterations (default = 100);}

\item{stand}{default standardize data TRUE}

\item{log}{default log-transformation TRUE}

\item{verbose}{default TRUE}

\item{fast}{default use fast RcppArmadillo implementation}
}
\value{
vector of performance measures with names indicating the model parameter
}
\description{
Perform the HCP normalization algorithm on a grid of model parameters
}
\details{
This function can be used to find the optimal model parameters
with a used-defined performance function
}
\examples{
\dontrun{
library(BiocParallel)
library(Rhcpp)
register(MulticoreParam(3))
kRange <- c(10, 20)
lambdaRange <- c(1, 5, 10, 20)
data(rhcppdata)
F <- rhcppdata$F
Y <- rhcppdata$Y
##not really meaning full performance function
res <- hcppcv(Z, Y, kRange, lambdaRange, performance=function(res) sum(res$Res))
hist(res)
which.min(res)
}
}
\author{
mvaniterson
}

