r_hcp <- function(Z, Y, k, lambda1, lambda2, lambda3, iter=100) {
    ## convergence criteria
    tol <- 1e-6

    A <- matrix(0, ncol(Z), k)
    W <- matrix(0, nrow(Z), k)
    n <- k*ncol(Y)
    ##B <- matrix(runif(n), ncol(Z), ncol(Y))
    B <- matrix((1:n)/n, k, ncol(Y))

    n1 <- nrow(Z)
    d1 <- ncol(Z)

    n2 <- nrow(Y)
    d2 <- ncol(Y)

    if(n1 != n2)
        message('number of rows in F and Y must agree')

    if (k < 1 | lambda1 < 0 | lambda2 < 0 | lambda3 < 0 )
        message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')


    ##predefine for slight preformance improvement
    diagB <- diag(k)
    diagW <- diag(k)
    diagA <- diag(nrow(A))
    U1 <- lambda1*solve(lambda1*crossprod(Z) + diagA*lambda3)%*%t(Z)

    if(iter > 0) {
        o <- numeric(iter)
        for(ii in 1:iter) {
            ##o[ii] <- norm(Y-W%*%B, type="F")^2 + lambda1*norm(W-Z%*%A, type="F")^2 + lambda2*norm(B, type="F")^2 + lambda3*norm(A, type="F")^2
            o[ii] <- sum((Y-W%*%B)^2) + sum((W-Z%*%A)^2)*lambda1 + sum(B^2)*lambda2 + lambda3*sum(A^2)           
            W <- (tcrossprod(Y, B) + lambda1*Z%*%A) %*% solve(tcrossprod(B) + lambda1*diagB)
            B <- solve(crossprod(W) + lambda2*diagW, crossprod(W,Y))            
            ##A <- solve(lambda1*crossprod(Z) + diagA*lambda3), lambda1*crossprod(Z,W))
            A <- U1%*%W
        
            if(ii > 1)  {
                if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
                    break
            }

        }
    }

    list(W=W, B=B, A=A, o=o, iter=ii)
}
