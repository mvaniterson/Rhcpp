##' Perform the HCP normalization algorithm on a grid of model parameters
##'
##' This function can be used to find the optimal model parameters
##' with a used-defined performance function
##'
##' @title Perform the HCP normalization algorithm on a grid of model parameters
##' @param Z a matrix nxd of known covariates, where n is the number of
##' subjects and d is the number of known covariates. *must be standardize
##' (columns have 0 mean and constant SS).
##' @param Y a matrix of nxg of expression data (must be standardized (columns
##' scaled to have constant SS and mean 0). ** use standardize function to standardize F and Y.
##' @param x vector of responses.
##' @param kRange multiple numbers of inferred hidden components (k is an integer)
##' @param lambdaRange multiple model parameters
##' @param iter (optional) iter: number of iterations (default = 100);
##' @param stand default standardize data TRUE
##' @param log default log-transformation TRUE
##' @param fast default use fast RcppArmadillo implementation
##' @param verbose default TRUE
##' @param performance function accepting res with res$Res the transformed Residuals
##' @return vector of performance measures with names indicating the model parameter
##' @importFrom BiocParallel bpworkers bplapply
##' @author mvaniterson
##' @examples
##' \dontrun{
##' library(BiocParallel)
##' library(Rhcpp)
##' register(MulticoreParam(3))
##' kRange <- c(10, 20)
##' lambdaRange <- c(1, 5, 10, 20)
##' data(rhcppdata)
##' F <- rhcppdata$F
##' Y <- rhcppdata$Y
##' ## we do not have response for this data
##' x <- rnorm(nrol(Y))
##' ##not really meaning full performance function
##' res <- hcppcv(F, Y, x, kRange, lambdaRange, performance=function(res) sum(res$Res))
##' res
##' }
hcppcv <- function(Z, Y, x, kRange=c(10, 20), lambdaRange=NULL, lambda1Range=NULL, lambda2Range=NULL, lambda3Range=NULL,
                   performance=NULL, iter=100, stand=TRUE, log=TRUE, verbose=TRUE, fast=TRUE) {

    if(is.null(performance))
        stop("A model performance function that accepts the output of hcp should be provided!")

    if(!is.null(lambda1Range) & !is.null(lambda2Range) & !is.null(lambda3Range))
        par <- expand.grid(k=kRange, lambda1=lambda1Range, lambda2=lambda2Range, lambda3=lambda3Range)
    else    
        par <- expand.grid(k=kRange, lambda1=lambdaRange, lambda2=lambdaRange, lambda3=lambdaRange)

    if(nrow(par) < bpworkers())
        stop("Number of workers:", bpworkers(), "should be smaller then the number of models to fit:", nrow(par))

    ##initial run perform log-transformation and standarization only once if necessary
    t0 <- proc.time()
    init <- hcpp(Z, Y, x, k = par$k[1], lambda1 =  par$lambda1[1], lambda2 = par$lambda2[1], lambda3 = par$lambda3[1], iter=iter, stand=stand, log=log, verbose=verbose, fast=fast)
    resinit <- performance(init)    
    estimatedTime <- (sum(par$k[-1])/par$k[1])*(50/init$niter)*(proc.time() - t0)[3]/bpworkers()

    if(verbose)
        message(paste0("Fitting all, ", nrow(par), ", models will approximately take: ", .sec2time(estimatedTime)))

    Y <- init$Y
    Z <- init$Z
    x <- (x - mean(x))/sqrt(sum((x - mean(x))^2))

    map <- function(i) {
        k <- par$k[i]
        lambda1 <- par$lambda1[i]
        lambda2 <- par$lambda2[i]
        lambda3 <- par$lambda3[i]
        message(paste("optimizing k = ", k, "lambda1 = ", lambda1, "lambda2 = ", lambda2, "lambda3 = ", lambda3))

        res <- hcpp(Z, Y, x, k = k, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, iter=iter, stand=FALSE, log=FALSE, verbose=verbose, fast=fast)
        performance(res)
    }

    gc() ##reduce memory footprint before running in parallel

    res <- bplapply(2:nrow(par), map)
    res <- c(resinit, unlist(res))
    names(res) <- apply(par, 1, paste0, collapse=":")
    res
}
