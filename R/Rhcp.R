##' pure R implementation of the HCP algorithm
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
##' @title pure R implementation of the HCP algorithm
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
##' @return list
##' Z: matrix of hidden components, dimensionality: nxk,
##' B: matrix of effects of hidden components, dimensionality: kxg,
##' o: value of objective function on consecutive iterations.
##' @author mvaniterson
##' @references \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0068141}
##' @examples
##' data(hcp)
##' F <- data$F
##' Y <- data$Y
##' k <- 10
##' lambda1 <- 20
##' lambda2 <- 1
##' lambda3 <- 1
##' iter <- 100
##' Rres <- r_hcp(Y, F, k, lambda1, lambda2, lambda3, iter)
r_hcp <- function(F, Y, k, lambda1, lambda2, lambda3, iter=NULL)
  {    
  
    if(is.null(iter))
      iter <- 100

    ## convergence criteria
    tol <- 1e-6

    U <- matrix(0, ncol(F), k)
    Z <- matrix(0, nrow(F), k)
    n <- ncol(Z)*ncol(Y)
    B <- matrix(runif(n), ncol(Z), ncol(Y))
    ## B <- matrix((1:n)/n, ncol(Z), ncol(Y))

    n1 <- nrow(F)
    d1 <- ncol(F)

    n2 <- nrow(Y)
    d2 <- ncol(Y)

    if(n1 != n2)
      message('number of rows in F and Y must agree')

    if (k < 1 | lambda1 < 1e-6 | lambda2 < 1e-6 | lambda3 < 1e-6 )
      message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')


    ##predefine for slight preformance improvement
    diagB <- diag(nrow(B))
    diagZ <- diag(ncol(Z))
    diagU <- diag(nrow(U))
    
    if(iter > 0)
      {
        o <- numeric(iter)
        for(ii in 1:iter)
          {
            ##o[ii] <- norm(Y-Z%*%B, type="F") + norm(Z-F%*%U, type="F")*lambda1 + norm(B, type="F")*lambda2 + lambda3*norm(U, type="F")
            o[ii] <- sum((Y-Z%*%B)^2) + sum((Z-F%*%U)^2)*lambda1 + sum(B^2)*lambda2 + lambda3*sum(U^2)

            ##Z <- (Y%*%t(B) + lambda1*F%*%U) %*% solve(B%*%t(B) + lambda1*diagB)
            Z <- (tcrossprod(Y, B) + lambda1*F%*%U) %*% solve(tcrossprod(B) + lambda1*diagB)
            
            ##B <- solve(t(Z)%*%Z + lambda2*diagZ, t(Z)%*%Y)
            B <- solve(crossprod(Z) + lambda2*diagZ, crossprod(Z,Y))

            ##U <- solve(t(F)%*%F*lambda1 + lambda3*diagU, lambda1*t(F)%*%Z)
            U <- solve(crossprod(F)*lambda1 + lambda3*diagU, lambda1*crossprod(F,Z))

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


    list(Z=Z, B=B, U=U, o=o)
  }


