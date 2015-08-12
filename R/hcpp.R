##' Rcpp implementation of the HCP algorithm
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
##' @param Z a matrix nxd of known covariates, where n is the number of
##' subjects and d is the number of known covariates. *must be standardize
##' (columns have 0 mean and constant SS).
##' @param Y a matrix of nxg of expression data (must be standardized (columns
##' scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
##' @param X vector of responses.
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
##' data(rhcppdata)
##' F <- rhcppdata$F
##' Y <- rhcppdata$Y
##' k <- 10
##' lambda1 <- 20
##' lambda2 <- 1
##' lambda3 <- 1
##' iter <- 100
##' ##and run
##' ##Rres <- hcpp(Y, F, k, lambda1, lambda2, lambda3, iter)
hcpp <- function(Z, Y, X, k, lambda1, lambda2, lambda3, iter=100, stand=TRUE, log=TRUE, fast=TRUE, verbose=TRUE) {
    t0 <- proc.time()
    if(is.null(Z)) {
        message("Assume only hidden components")
        Z <- matrix(runif(2*nrow(Y), 1, 10), nrow=nrow(Y), ncol=2)
        lambda3 <- 0 ##do not penalized the coefficients with an effect on the known covariates
    }
    ## (0) check dimensions
    if(nrow(Y) != nrow(Z) | nrow(Y) != nrow(X))
        stop("Rows represent the samples for both Y, Z and X!")

    if(sum(is.na(Z)) > 0 | sum(is.na(Y)) > 0 | sum(is.na(X)) > 0)
        stop("NA's in input data are not allowed!")
 
    ## (1) take log of read counts
    if(log) {
        message("Log-transformation data...")
        Y <- log2(2+Y)
        if(sum(is.na(Y)) > 0)
            stop("NA's introduced by log-transformation these are not allowed!")
    }

    standardize <- function(x) {
        x <- x - outer(rep(1, nrow(x)), colMeans(x)) ##center
        x <- x*outer(rep(1, nrow(x)), 1/sqrt(colSums(x^2))) ##scale SS
        ##x <- apply(x, 2, function(x) x - mean(x))
        ##x <- apply(x, 2, function(x) x/sqrt(sum(x^2)))
        x
    }

    ## (2) standardize the data
    if(stand) {
        message("Standardize data...")
        Y <- standardize(Y)
        Z <- standardize(Z)
        X <- standardize(X)
    }

    if(sum(is.na(Z)) > 0 | sum(is.na(Y)) > 0 | sum(is.na(X)) > 0) {
        message(paste("Row(s) (Z):", paste(which(apply(Z, 1, function(x) any(is.na(x)))), collapse=", ")))
        message(paste("Row(s) (Y):", paste(which(apply(Y, 1, function(x) any(is.na(x)))), collapse=", ")))
        message(paste("Row(s) (X):", paste(which(apply(X, 1, function(x) any(is.na(x)))), collapse=", ")))
        stop("NA's introduced by the standardization these are not allowed!")
    }

    ## (3) HCP
    if(fast) {
        message("Run RcppArmadillo implemented HCP algorithm...")
        res <- rcpparma_hcpp(Z, Y, X, k, lambda1, lambda2, lambda3, iter)
    }
    else {
        message("Run plain R implemented HCP algorithm...")
        res <- r_hcpp(Z, Y, X, k, lambda1, lambda2, lambda3, iter)
    }

    if(verbose)
        message(paste("Finished after", res$niter, "iterations of", iter, "iterations."))

    W <- res$W
    B <- res$B
    o <- res$o
    G <- res$G
    g <- as.vector(G[,1])
    niter <- res$niter
    err <- as.vector(sqrt(colSums((Y - X%*%G - W%*%B)^2)/(nrow(Y)- ncol(W) - ncol(X) - 1)))
    ##by default inference on first column of X
    pval <- 2*pnorm(-abs(g/err)) ##approximation to t; n is usually large enough
    names(pval) <- colnames(Y)

    if(verbose)
        message(paste("The batch correction took:", round((proc.time() - t0)[3], 2), "seconds."))

    return(list(Res = Y - X%*%G - W%*%B, Cov=W, B=B, o=o, gamma=g, err=as.vector(err), pval=as.vector(pval), Y=Y, Z=Z, X=X, niter=niter))
}
