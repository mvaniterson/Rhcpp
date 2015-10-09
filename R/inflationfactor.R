##' Gibbs sampler that fits normal mixtures
##' details follow
##' TODO improve implementation e.g. S3 to allow easy plotting etc.
##' @title Gibbs sampler
##' @param y numeric vector containing the observed mixture values
##' @param k number of components of the mixture default is 3 
##' @param niter length of the Gibbs sampling run default 5000
##' @param burnin length of the Burn-in default is 2000
##' @param verbose print progress information default is TRUE
##' @param plot plot traces progress default FALSE
##' @param alpha priors for variances inverse-Gamma distribution
##' @param beta priors for variances inverse-Gamma distribution 
##' @param lam priors for means Normal distribution
##' @param tau priors for means Normal distribution
##' @param g priors for proportions Dirichlet distribution
##' @param rho1 prior for variances Beta distribution
##' @param rho2 prior for variances Beta distribution
##' @return posterior samples from the fitted mixture distribution
##' @author mvaniterson
##' @references Implementation is based on a version from Zhihui Liu \url{https://macsphere.mcmaster.ca/handle/11375/9368} 
##' @export
##' @examples
##' theta <- c(0.9, 1, 0, 3, 1)
##' y <- rnormmix(5000, theta)
##' gs <- gibbssampler(y, niter=500, burnin=200)
##' plotTraces(gs)
##' plotCI(gs)
##' points(1, 0.9, pch=16, col=2, cex=2)
##' theta.hat <- apply(gs, 2, median)
##' plotnormmix(y, theta.hat, n=100)
gibbssampler <- function(y, k=3, niter=5000, burnin=2000, verbose=TRUE, plot=FALSE,
                         alpha=1.28, beta=0.36*var(y), 
                         lam = c(0, 1, -1), tau =  c(1000, rep(100, 2)), 
                         g = c(90,5,5), 
                         rho1=0.1842, rho2=3.5) {

    n <- length(y)
    mu <- lam           ##sample(y, k)
    sig <- rep(sd(y)/k, k)
    p <- g/sum(g)       ##rep(1/k, k)
    nj <- sj <- sj2 <- rep(0, k)
    mixparam <- list(p = p, mu = mu, sig = sig)
    gibbsmu <- gibbssig <- gibbsp <- matrix(0, nrow= niter, ncol = k)
    for(i in 1:niter) {

        prob <- apply(matrix(y, ncol=1), 1, function(x) mixparam$p*dnorm(x, mean = mixparam$mu, sd = mixparam$sig))
        z <- apply(prob, 2, function(x) sample(1:k, 1, prob=x))

        for(j in 1:k) {
            nj[j] <- sum(z == j)
            sj[j] <- sum(as.numeric(z == j)*y)
        }

        ## repeat {
        ##     gibbsmu[i,] <- rnorm(k, mean = (lam*tau + sj)/(nj + tau), sd = sqrt(mixparam$sig^2/(tau+nj)))
        ##     if(max(gibbsmu[i,]) < max(y) & min(gibbsmu[i,]) > min(y))
        ##         break
        ## }

        gibbsmu[i,] <- rnorm(k, mean = (lam*tau + sj)/(nj + tau), sd = sqrt(mixparam$sig^2/(tau+nj)))

        mixparam$mu <- gibbsmu[i,]

        for(j in 1:k) {
            sj2[j] <- sum(as.numeric(z == j) * (y - mixparam$mu[j])^2)
        }

        gibbssig[i, ] <- sqrt(.rigamma(k, alpha + 0.5*(nj + 1), beta + 0.5*tau *(mixparam$mu -lam)^2 + 0.5*sj2))
        mixparam$sig <- gibbssig[i,]

        if(k == 2)
            gibbsp[i,] <- c(rbeta(1, rho1+nj[1], rho2+nj[2]), 1 - rbeta(1, rho1+nj[1], rho2+nj[2])) ##informative beta prior
        else
            gibbsp[i,] <- .rdirichlet(1, par = nj + g)
        mixparam$p <- gibbsp[i,]

        if(verbose & (i %% as.integer(niter/10) == 0))
            message(paste("At iteration", i, "of", niter, "..."))

        if(plot & (i %% as.integer(niter/10) == 0)) {
            gs <- data.frame(p = gibbsp, mu = gibbsmu, sigma = gibbssig)
            plotTraces(gs[1:i,])
        }

    }

    gs <- data.frame(p = gibbsp, mu = gibbsmu, sigma = gibbssig)

    if(!is.null(burnin))
        return(gs[seq(burnin, niter),])
    else
        return(gs)
}

##' Show traces of the gibbs sampling
##'
##' details follow not really generic
##' @title plot traces of the samples drawn from the conditional distribtions
##' @param gs gibssampler object
##' @param theta optionally select parameters to plot
##' @return trace plot of the different parameters
##' @author mvaniterson
plotTraces <- function(gs, theta=rep(NA, ncol(gs))) {    
    op <- par(mfcol=c(ncol(gs)/3, 3), mar=c(2,4,2,2))
    for(i in 1:ncol(gs)) {
        plot(gs[,i], ylab=colnames(gs)[i], type="l", xlab="", main = "", lwd=0.3)
        abline(h=theta[i], col=2, lwd=2)
    }
    par(op)
}

##' plot normal based confidence ellipse on two parameters
##'
##' details follow not really generic
##' @title plot confidence ellipse
##' @param gs gibssampler object
##' @param theta select parameters to plot
##' @param alphas alfa level(s) for the confidence ellipses default 0.95, 0.9 and 0.75.
##' @param ... arguments to ellipse
##' @return plot
##' @author mvaniterson
##' @importFrom ellipse ellipse
plotCI <- function(gs, theta = c("sigma.1", "p.1"), alphas=c(0.95, 0.9, 0.75), ...) {
    d <- gs[,theta]
    plot(d, pch=20, xlab=expression(lambda), ylab=expression(1-epsilon), bty='n',
         main=c("median at:", round(c(median(gs[,theta[1]]), median(gs[,theta[2]])),3)))
    points(median(gs[,theta[1]]), median(gs[,theta[2]]), col=3, pch=17, cex=2)
    for(alpha in alphas)
        lines(ellipse(cov(d), centre=colMeans(d), level=alpha), col="blue", ...)
}

##' sample from a normal mixture
##'
##' details follow
##' @title sample from a normal mixture
##' @param n size 
##' @param theta parameters
##' @return n samples from a normal mixture with parameters theta
##' @author mvaniterson
##' @export
rnormmix <- function(n, theta){
    epsilon1 <- theta[1]
    lambda <- theta[2]
    mu1 <- theta[3]
    sigma1 <- theta[4]
    sigma2 <- theta[5]

    U <- runif(n)
    n1 <- sum(U < epsilon1)
    n2 <- sum(U >=  epsilon1)
    mu <- rnorm(n2, mu1, sigma1)
    c(rnorm(n1, 0, lambda), mu + rnorm(n2, 0, sigma2))
}

##' density of a k-component normal mixture
##'
##' details follow
##' @title density of a k-component normal mixture
##' @param x x like dnorm(x, ...
##' @param theta parameters of the mixture proportion, mean and sd
##' @return density of a k-component normal mixture
##' @author mvaniterson
##' @export
dnormmix <- function(x, theta){
    if(length(theta) %% 3 != 0)
        stop("Length of theta should be a multiple of three!")
    ncomp <- length(theta)/3
    y <- 0*x
    for(k in 1:ncomp)
        y <- y + theta[k]*dnorm(x, mean=theta[k+3], sd=theta[k+6])
    y
}

##' plot normal mixtures
##'
##' details follow
##' @title plot normal mixtures
##' @param data vector of test statistics
##' @param theta parameters describing the mixture components
##' @param ... arguments passed to hist
##' @return return plot with histogram of the data and mixture and individual components
##' @author mvaniterson
##' @export
plotnormmix <- function(data, theta, ...) {
    if(length(theta) %% 3 != 0)
        stop("Length of theta should be a multiple of three!")
    hist(data, freq=FALSE, ...)
    f <- function(x) dnormmix(x, theta)
    curve(f, add=TRUE, col=1, lwd=2)
    ncomp <- length(theta)/3
    for(k in 1:ncomp) {
        f <- function(x) theta[k]*dnorm(x, mean=theta[k+3], sd=theta[k+6])
        curve(f, add=TRUE, col=k+1, lwd=2)
    }
}
