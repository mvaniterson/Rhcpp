r_hcpp1 <- function(Z, Y, x, k, lambda, iter=500, verbose, tol = 1e-6) {

    if(nrow(Z) != nrow(Y))
        message('number of rows in Z and Y must agree')

    if (k < 1 | lambda < 0)
        message('lambda1 must be positive and/or k must be an integer')

    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(Z)

    W <- matrix(runif(n*k), n, k)
    A <- matrix(runif(q*k), q, k)
    B <- matrix(runif(k*p), k, p)
    G <- matrix(runif(1*p), 1, p)

    o <- numeric(iter)
    for(ii in 1:iter) {
        o[ii] <- sum((Y - x%*%G - W%*%B)^2) + lambda*sum((W - Z%*%A)^2) + lambda*(sum(B^2)) ##frobenius norm
        A <- qr.solve(crossprod(Z), crossprod(Z, W))
        B <-  qr.solve(crossprod(W) + lambda*diag(k), crossprod(W, Y - x%*%G))
        G <- qr.solve(crossprod(x), crossprod(x, Y - W%*%B))
        WT <- qr.solve(t(tcrossprod(B) + lambda*diag(k)), t(tcrossprod(Y - x%*%G, B) + lambda*Z%*%A))
        W <- t(WT)

        if(ii > 1) {
            if(abs(o[ii] - o[ii-1])/o[ii] < tol)
                break
        }

        if(ii %% 20 == 0 & verbose)
            message(paste("at iteration:", ii, "of max", iter, "(rel. diff.", signif(abs(o[ii] - o[ii-1])/o[ii], 4), ")"))
    }

    list(W=W, B=B, A=A, o=o, G=G, niter=ii)
}

r_hcpp3 <- function(Z, Y, X, k, lambda1, lambda2, lambda3, iter=100, verbose, tol=1e-6) {

    W <- matrix(0, nrow(Z), k)
    A <- matrix(0, ncol(Z), k)
    G <- matrix(0, ncol(X), ncol(Y))

    n <- k*ncol(Y)
    B <- matrix((1:n)/n, k, ncol(Y)) ##B <- matrix(runif(n), ncol(Z), ncol(Y))

    if(nrow(Z) != nrow(Y))
        message('number of rows in Z and Y must agree')

    if (k < 1 | lambda1 < 0 | lambda2 < 0 | lambda3 < 0 )
        message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')

    ##predefine
    diagB <- diagW <- diag(k)
    diagA <- diag(nrow(A))

    o <- numeric(iter)
    for(ii in 1:iter) {
        o[ii] <- norm(Y-X%*%G-W%*%B, type="F")^2 + lambda1*norm(W-Z%*%A, type="F")^2 + lambda2*norm(B, type="F")^2 + lambda3*norm(A, type="F")^2
        W <- (tcrossprod(Y-X%*%G, B) + lambda1*Z%*%A)%*%solve(tcrossprod(B) + lambda1*diagB)
        B <- solve(crossprod(W) + lambda2*diagW, crossprod(W,Y-X%*%G))
        A <- solve(crossprod(Z) + (lambda3/lambda1)*diagA, crossprod(Z,W))
        G <- solve(crossprod(X), crossprod(X, Y-W%*%B))

        if(ii > 1) {
            if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
                break
        }
    }

    list(W=W, B=B, A=A, o=o, G=G, niter=ii)
}
