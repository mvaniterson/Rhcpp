r_hcp <- function(F, Y, k, lambda1, lambda2, lambda3, iter=NULL)
  {

    if(is.null(iter))
      iter <- 100

    ## convergence criteria
    tol <- 1e-6

    U <- matrix(0, ncol(F), k)
    Z <- matrix(0, nrow(F), k)
    n <- k*ncol(Y)
    ## B <- matrix(runif(n), ncol(Z), ncol(Y))
    B <- matrix((1:n)/n, k, ncol(Y))

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


