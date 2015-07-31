r_hcpp <- function(Z, Y, x, k, lambda1, lambda2, lambda3, iter=100) {

    ## convergence criteria
    tol <- 1e-6

    W <- matrix(0, nrow(Z), k)
    A <- matrix(0, ncol(Z), k)    
    g <- matrix(0, 1, ncol(Y))
    
    n <- k*ncol(Y)    
    B <- matrix((1:n)/n, k, ncol(Y)) ##B <- matrix(runif(n), ncol(Z), ncol(Y))

    if(nrow(Z) != nrow(Y))
        message('number of rows in F and Y must agree')

    if (k < 1 | lambda1 < 1e-6 | lambda2 < 1e-6 | lambda3 < 1e-6 )
        message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')

    ##predefine 
    diagB <- diagW <- diag(k)    
    diagA <- diag(nrow(A))
    xtx <- as.numeric(crossprod(x)) ##since this is a scalar
    
    if(iter > 0) {
        o <- numeric(iter)
        for(ii in 1:iter) {
            o[ii] <- norm(Y-W%*%B, type="F")^2 + lambda1*norm(W-Z%*%A, type="F")^2 + lambda2*norm(B, type="F")^2 + lambda3*norm(A, type="F")^2            
            W <- (tcrossprod(Y-x%*%gamma, B) + lambda1*Z%*%A)%*%solve(tcrossprod(B) + lambda1*diagB)
            B <- solve(crossprod(W) + lambda2*diagW, crossprod(W,Y-x%*%gamma))         
            A <- solve(crossprod(F) + (lambda3/lambda1)*diagA, crossprod(F,Z))                       
            g <- t(x)%*%(Y - Z%*%B)/xtx
            
            if(ii > 1) {
                if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
                    break
            }
        }
    }
    
    list(W=W, B=B, A=A, o=o, g=g, iter=ii)
}
