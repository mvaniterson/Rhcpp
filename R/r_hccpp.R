r_hcpp <- function(F, Y, x, k, lambda1, lambda2, lambda3, iter=NULL)
  {

    if(is.null(iter))
      iter <- 100

    ## convergence criteria
    tol <- 1e-6

    U <- matrix(0, ncol(F), k)
    Z <- matrix(0, nrow(F), k)
    n <- k*ncol(Y)
    ##B <- matrix(runif(n), ncol(Z), ncol(Y))
    B <- matrix((1:n)/n, k, ncol(Y))
    gamma <- matrix(0, 1, ncol(Y))

    n1 <- nrow(F)
    d1 <- ncol(F)

    n2 <- nrow(Y)
    d2 <- ncol(Y)

    if(n1 != n2)
      message('number of rows in F and Y must agree')

    if (k < 1 | lambda1 < 1e-6 | lambda2 < 1e-6 | lambda3 < 1e-6 )
      message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')


    ##predefine for slight preformance improvement
    diagB <- diag(k)
    diagZ <- diag(k)
    diagU <- diag(nrow(U))

    if(iter > 0)
      {
        o <- numeric(iter)
        for(ii in 1:iter)
          {
            ##o[ii] <- norm(Y-Z%*%B, type="F") + norm(Z-F%*%U, type="F")*lambda1 + norm(B, type="F")*lambda2 + lambda3*norm(U, type="F")
            o[ii] <- sum((Y-x%*%gamma-Z%*%B)^2) + sum((Z-F%*%U)^2)*lambda1 + sum(B^2)*lambda2 + lambda3*sum(U^2)

            ##Z <- (Y%*%t(B) + lambda1*F%*%U) %*% solve(B%*%t(B) + lambda1*diagB)
            Z <- (tcrossprod(Y-x%*%gamma, B) + sqrt(lambda1)*F%*%U) %*% solve(tcrossprod(B) + lambda1*diagB)

            ##B <- solve(t(Z)%*%Z + lambda2*diagZ, t(Z)%*%Y)
            B <- solve(crossprod(Z) + lambda2*diagZ, crossprod(Z,Y-x%*%gamma))
            ##can this be simplified??
            ##B <- solve(crossprod(Z) + lambda2*diagZ)%*%crossprod(Z,Y)

            ##U <- solve(t(F)%*%F*lambda1 + lambda3*diagU, lambda1*t(F)%*%Z)
            ##U <- solve(crossprod(F)*lambda1 + lambda3*diagU, lambda1*crossprod(F,Z))
            U <- solve(crossprod(F) + diagU*(lambda3/lambda1), crossprod(F,Z))

            gamma <- solve(crossprod(x))%*%crossprod(x, Y - Z%*%B)
                
            if(ii > 1)
              {
                if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
                  break
              }

          }
      }

    ##erroro <- sum(sum((Y-Z*B).^2))./sum(sum(Y.^2)) + sum(sum((Z-F*U).^2))./sum(sum((F*U).^2))
    ##error1 <- sum(sum((Y-Z*B).^2))./sum(sum(Y.^2))
    ##error2 <- sum(sum((Z-F*U).^2))./sum(sum((F*U).^2))

    ##dz <- Z*(B*t(B) + lambda1*eye(size(B,1)))-(Y*t(B) + lambda1*F*U)
    ##db <- (t(Z)*Z + lambda2*eye(size(Z,2)))*B - t(Z)*Y
    ##du <- (t(F)*F*lambda1 + lambda3*eye(size(U,1)))*U-lambda1*t(F)*Z

    list(Z=Z, B=B, U=U, o=o, gamma=gamma, iter=ii)
  }

hcpp <- function(F, Y, x, k, lambda1, lambda2, lambda3, iter=NULL, stand=TRUE, log=TRUE, fast=TRUE, verbose=TRUE)
  {
    t0 <- proc.time()
    if(is.null(F)) {
        message("Assume only hidden components")
        F <- matrix(runif(2*nrow(Y), 1, 10), nrow=nrow(Y), ncol=2)
        lambda3 <- 0 ##do not penalized the coefficients with an effect on the known covariates
    }
    ## (0) check dimensions
    if(nrow(Y) != nrow(F))
      stop("Rows represent the samples for both Y and F!")

    if(sum(is.na(F)) > 0 | sum(is.na(Y)) > 0)
      stop("NA's in input data are not allowed!")

    ## (1) take log of read counts
    if(log)
      {
        message("Log-transformation data...")
        Y <- log2(2+Y)
        if(sum(is.na(Y)) > 0)
          stop("NA's introduced by log-transformation these are not allowed!")
      }

    standardize <- function(x)
      {
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
      F <- standardize(F)
      x <- (x - mean(x))/sqrt(sum((x - mean(x))^2))
    }

    if(sum(is.na(F)) > 0 | sum(is.na(Y)) > 0) {        
      message(paste("Row(s) (F):", paste(which(apply(F, 1, function(x) any(is.na(x)))), collapse=", ")))
      message(paste("Row(s) (Y):", paste(which(apply(Y, 1, function(x) any(is.na(x)))), collapse=", ")))
      stop("NA's introduced by the standardization these are not allowed!")
     }

    ## (3) HCP
    if(fast) {
      message("Run RcppArmadillo implemented HCP algorithm...")
      res <- rcpparma_hcp(F, Y, k, lambda1, lambda2, lambda3, iter)
    }
    else
      {
        message("Run plain R implemented HCP algorithm...")
        res <- r_hcpp(F, Y, x, k, lambda1, lambda2, lambda3, iter)
      }

    if(verbose)
      message(paste("Finished after", res$iter, "iterations of", iter, "iterations."))

    Z <- res$Z
    B <- res$B
    gamma <- res$gamma
    ##should I force x to have dim nx1?
    err <- as.vector(sqrt(colSums((Y - Z%*%B - x%*%gamma)^2)/(nrow(Y)-2)))
    pval <- 2*pnorm(-abs(gamma/err)) ##approximation to t; n is usually large enough
    names(pval) <- colnames(Y)
    
    if(verbose)
        message(paste("The batch correction took:", round((proc.time() - t0)[3], 2), "seconds."))
    
    return(list(Res = Y - Z%*%B, Cov=Z, B=B, gamma=gamma, err=err, pval=pval, Y=Y, F=F))
  }

hcppcv <- function(F, Y, x, kRange=c(10, 20), lambdaRange=c(1, 5, 10, 20), performance=NULL, iter=100, stand=TRUE, log=TRUE, verbose=TRUE, fast=TRUE) {

    if(is.null(performance))
        stop("A model performance function that accepts the output of hcp should be provided!")

    par <- expand.grid(k=kRange, lambda1=lambdaRange, lambda2=lambdaRange, lambda3=lambdaRange)

    ##initial run perform log-transformation and standarization only once if necessary
    t0 <- proc.time()
    init <- hcpp(F, Y, x, k = par$k[1], lambda1 =  par$lambda1[1], lambda2 = par$lambda2[1], lambda3 = par$lambda3[1], iter=iter, stand=stand, log=log, verbose=verbose, fast=fast)
    resinit <- performance(init)
    estimatedTime <- nrow(par)*(proc.time() - t0)[3]/bpworkers()
    
    if(verbose) 
        message(paste0("Fitting all, ", nrow(par), ", models will approximately take: ", .sec2time(estimatedTime)))

    Y <- init$Y
    F <- init$F

    map <- function(i) {
        k <- par$k[i]
        lambda1 <- par$lambda1[i]
        lambda2 <- par$lambda2[i]
        lambda3 <- par$lambda3[i]
        message(paste("optimizing k =", k, "lambda1 = ", lambda1, "lambda2 = ", lambda2, "lambda3 = ", lambda3))

        res <- hcpp(F, Y, x, k = k, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, iter=iter, stand=FALSE, log=FALSE, verbose=verbose, fast=fast)
        performance(res)
    }
    
    gc() ##reduce memory footprint before running in parallel
    
    res <- bplapply(2:nrow(par), map)
    res <- c(resinit, unlist(res))
    names(res) <- apply(par, 1, paste0, collapse=":")
    res
}

