##' R/Rcpp implementation of the HCP algorithm
##'
##' Objective:
##' This function solves the following problem:
##' argmin_{Z,B,U}   ||Y-Z*B||_2 + lambda1*||Z-F*U||_2 + lambda2*||B||_2 +
##' lambda_3||U||_2
##'
##' To use the residual data: Residual = Y - Z*B
##'
##' Note: k>0, lambda1>0, lambda2>0, lambda3>0 must be set by the user based on
##' the data at hand. one can set these values using cross-validation, by
##' evaluating the "performance" of the  resulting residual data on a desired
##' task. typically, if lambda>5, then hidden factors match the known covariates closely.
##' @title R/Rcpp implementation of the HCP algorithm
##' @param F a matrix nxd of known covariates, where n is the number of
##' subjects and d is the number of known covariates. *must be standardize
##' (columns have 0 mean and constant SS).
##' @param Y a matrix of nxg of expression data (must be standardized (columns
##' scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
##' @param k  number of inferred hidden components (k is an integer)
##' @param lambda1 model parameter 1
##' @param lambda2 model parameter 2
##' @param lambda3 model parameter 3
##' @param iter (optional) iter: number of iterations (default = 100);
##' @param stand default standardize data TRUE
##' @param log default log-transformation TRUE
##' @param fast default use fast RcppArmadillo implementation
##' @param verbose default TRUE
##' @return list
##' Z: matrix of hidden components, dimensionality: nxk,
##' B: matrix of effects of hidden components, dimensionality: kxg,
##' o: value of objective function on consecutive iterations.
##' @author mvaniterson
##' @references \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0068141}
##' @examples
##' data(hcpp)
##' F <- hcpp$F
##' Y <- hcpp$Y
##' k <- 10
##' lambda1 <- 20
##' lambda2 <- 1
##' lambda3 <- 1
##' iter <- 100
##' ##and run
##' ##Rres <- hcp(Y, F, k, lambda1, lambda2, lambda3, iter)
hcp <- function(F, Y, k, lambda1, lambda2, lambda3, iter=NULL, stand=TRUE, log=TRUE, fast=TRUE, verbose=TRUE)
  {
    t0 <- proc.time()

    ## (0) check dimensions
    if(nrow(Y) != nrow(F))      
      stop("Rows represent the samples for both Y and F!")
      
    ## (1) take log of read counts
    if(log)
      {
        message("Log-transformation data...")
        Y <- log2(2+Y)
      }

    if(sum(is.na(F)) > 0 | sum(is.na(Y)) > 0)
      stop("1: NA's in input data are not allowed!")

    standardize <- function(x)
      {
        x <- apply(x, 1, function(x) x - mean(x))
        x <- apply(x, 1, function(x) x/sqrt(sum(x^2)))
        x
      }

    ## (2) standardize the data
    if(stand) {
      message("Standardize data...")
      Yn <- scale(Y, scale=apply(Y, 2, function(x) sqrt(sum((x - mean(x))^2))))
      Fn <- scale(F, scale=apply(F, 2, function(x) sqrt(sum((x - mean(x))^2))))
      attributes(Yn)$`scaled:center` <- attributes(Yn)$`scaled:scale` <- NULL
      attributes(Fn)$`scaled:center` <- attributes(Fn)$`scaled:scale` <- NULL
    }
    else
      {
        Yn <- Y
        Fn <- F
      }

    if(sum(is.na(Fn)) > 0 | sum(is.na(Yn)) > 0)
      stop("2: NA's in input data are not allowed!")

    ## (3) HCP
    if(fast) {
      message("Run RcppArmadillo implemented HCP algorithm...")
      res <- rcpparma_hcp(Fn, Yn, k, lambda1, lambda2, lambda3, iter)      
    }
    else
      {
        message("Run plain R implemented HCP algorithm...")
        res <- r_hcp(Fn, Yn, k, lambda1, lambda2, lambda3, iter)
      }

    if(verbose)
      message(paste("Finished after", res$iter, "iterations of", iter, "iterations."))

    Z <- res$Z
    B <- res$B

    message(paste("The batch correction took:", round((proc.time() - t0)[3], 2), "seconds."))
    return(list(Res = t(Yn - Z%*%B), Cov = Z))
  }
