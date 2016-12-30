##SOME INITIAL INVESTIGATION OF THE HEFT METHOD
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3904522/

if(FALSE) {

    standardize <- function(x) {
        x <- x - outer(rep(1, nrow(x)), colMeans(x)) ##center
        x <- x*outer(rep(1, nrow(x)), sqrt(nrow(x)-1)/sqrt(colSums(x^2))) ##scale SS
        ##x <- apply(x, 2, function(x) x - mean(x))
        ##x <- apply(x, 2, function(x) x/sqrt(sum(x^2)))
        x
    }


    ##EM:
    Y <- Y
    Z <- Z
    X <- x

    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(Z)
    k <- 5

    alpha <- matrix(0, q, p)
    F <- matrix(0, k, p)

    Psi <- diag(n)
    Omega <- matrix(0, n, 1+q+k)
    Gamma <- matrix(0, 1+q+k, p)

    ## 1) initialize
    diag(Psi) <- 0.5
    beta <- matrix(0, 1, p)
    Lambda <- matrix(rnorm(n*k, 0, 1), n, k)

    Gamma[1,] <- beta
    Gamma[2:(1+q), ] <- alpha
    Gamma[(q+2):(1+q+k),] <- F

    iPsi <- 2*Psi

    ## 2/3) standardize
    H <- standardize(Y)
    Z <- standardize(Z)
    X <- standardize(X)

    Omega[,1] <- X
    Omega[,2:(1+q)] <- Z
    Omega[,(q+2):(1+q+k)] <- Lambda


    trHHT <- sum(diag(tcrossprod(H)))
    trHHT

    niter <- 100
    loglik <- numeric(niter)
    for(i in 1:niter) {

        ## 4) conditional variance Gamma given Y
        V <- solve(diag(1+q+k) + t(Omega)%*%iPsi%*%Omega)

        ## 5) conditional expectation Gamma given Y
        E <- V%*%t(Omega)%*%iPsi%*%H

        ## 6) Set the corresponding row of E to beta/alpha
        alpha <- solve(t(Z)%*%iPsi%*%Z, t(Z)%*%iPsi%*%(H - X%*%beta - Lambda%*%F))
        beta <- solve(t(X)%*%iPsi%*%X + 1, t(X)%*%iPsi%*%(H - Z%*%alpha - Lambda%*%F))
        E[1,] <- beta
        E[2:(1+q), ] <- alpha

        ## 7) Keep the fixed effect and the known covariates in the Omega matrix fixed, update the rest using:

        EGamma <- tcrossprod(E) - p*V
        T <- H%*%t(E)%*%solve(EGamma)
        Omega[,(q+2):(1+q+k)] <- T[,(q+2):(1+q+k)]

        ## 8/9)
        iPsi <- diag(n)/(trHHT/(n*p) - sum(diag(H%*%t(E)%*%t(Omega)))/(n*p))
        print(iPsi[1,1])

        ##update
        F <- E[(q+2):(1+q+k),]
        Lambda <- Omega[,(q+2):(1+q+k)]

        Gamma[1,] <- beta
        Gamma[2:(1+q), ] <- alpha
        Gamma[(q+2):(1+q+k),] <- F

        loglik[i] <- -0.5*iPsi[1,1]*sum(diag(tcrossprod(Y-Omega%*%Gamma))) - crossprod(beta) - sum(diag(tcrossprod(F))) - n*p*log(iPsi[1,1])

        loglik[i] <- -0.5*iPsi[1,1]* (trHHT - 2*sum(diag(tcrossprod(Omega%*%Gamma, H))) + 2*sum(diag(tcrossprod(Omega%*%Gamma)))) - tcrossprod(beta) - sum(diag(tcrossprod(F))) - n*p*log(iPsi[1,1])

        print(loglik[i])

        if(i > 1)
            if(abs(loglik[i]-loglik[i-1])/loglik[i] < 0.001)
                break
    }


    ##maximum likelihood
    ##HEFT
    Y <- t(counts[,id])
    Z <- covariates[id,-c(6, 9)]
    X <- matrix(response[id], ncol=1)

    n <- nrow(Y)
    p <- ncol(Y)
    q <- ncol(Z)
    k <- 5

    ## 1) standardize
    H <- standardize(Y)
    Z <- standardize(Z)
    X <- standardize(X)

    ## 2) initialize
    Omega <- matrix(0, n, 1+q+k)
    Gamma <- matrix(0, 1+q+k, p)
    Lambda <- matrix(rnorm(n*(1+q+k), 0, 1), n, 1+q+k)

    Omega[,1] <- X
    Omega[,2:(1+q)] <- Z

    Gamma[1, ] <- solve(crossprod(X), crossprod(X, H)) ##beta
    Gamma[2:(q+1),] <- solve(crossprod(Z), crossprod(Z, H- X%*%Gamma[1,]))  ##alpha
    Gamma[(q+2):(1+q+k),] <- matrix(0, k, p) ##F

    ## 3) iterate
    niter <- 25
    loglik <- numeric(niter)
    for(i in 1:niter) {
        Omega[,(q+2):(1+q+k)] <- Lambda[,(q+2):(1+q+k)]
        Gamma <- solve(crossprod(Omega) + diag(1+q+k), crossprod(Omega, H)) #ridge estimator for beta and F
        beta <- Gamma[1,,drop=FALSE]
        F <- Gamma[(q+2):(1+q+k),]
        Gamma[2:(q+1),] <- solve(crossprod(Z), crossprod(Z, H - Omega[,1]%*%beta - Omega[,(q+2):(1+q+k)]%*%F)) ##replace alpha with fixed estimator
        Lambda <- t(solve(tcrossprod(Gamma), tcrossprod(Gamma, H)))
        sigma2 <- sum((H - Omega%*%Gamma)^2)/(n*p) ##diag(sum(crossprod(A))) = sum(A^2)
        loglik[i] <- - 0.5*sum((H - Omega%*%Gamma)^2) - sum(beta^2) - sum(F^2)
        print(sigma2)
        print(loglik[i])
        if(i > 1){
            if(loglik[i] < loglik[i-1])
                break
        }
    }


    heft <- function(Y, Z, X, k, lambda=1, niter=50, tolerance=0.01) {
        n <- nrow(Y)
        p <- ncol(Y)
        q <- ncol(Z)
        m <- 1+q+k

        H <- standardize(Y)
        Z <- standardize(Z)
        X <- standardize(X)

        F <- matrix(0, k, p)
        Omega <- matrix(0, n, m)
        EGamma <- matrix(0, m, p)

        Omega[,1] <- X
        Omega[,2:(1+q)] <- Z
        Omega[,(q+2):m] <- matrix(rnorm(n*k, 0, 1), n, k)
        sigma2 <- 0.5

        HHT <- tcrossprod(H)

        loglik <- numeric(niter)
        for(i in 1:niter) {

            VG <- solve(lambda*diag(m) + crossprod(Omega, iPsi%*%Omega))
            EG <- VG%*%crossprod(Omega, iPsi%*%H)

            EG[1,] <- crossprod(X, (H - Omega[,-1]%*%EGamma[-1,]))/as.numeric((crossprod(X)+lambda))
            EG[2:(q+1),] <- solve(crossprod(Z), crossprod(Z, (H - Omega[,-c(2:(q+1))]%*%EGamma[-c(2:(q+1)),])))

            EGG <- tcrossprod(EG)+p*VG
            Omega <- H%*%t(EG)%*%solve(EGG)
            Omega[,1] <- X
            Omega[,2:(1+q)] <- Z

            S <- HHT/p - Omega%*%EG%*%t(H)/p
            sigma2 <- sum(diag(S))/n
            iPsi <- diag(n)/sigma2

            loglik[i] <- 0.5*log(det(iPsi)) - 0.5*sum(diag(tcrossprod(H - Omega%*%EG)%*%iPsi)) - sum(EG[1,]^2) - sum(EG[(q+2):m,]^2) ##-0.5*n*p*log(2*pi)

            print(paste0("iteration:", i))
            print(sigma2)
            print(loglik[i])

            ##mininimal number of iterations
            ##stop if difference is smaller or start to decreasing
            if(i > 10) {
                if(loglik[i] - loglik[i-1] < tolerance)
                    break
            }
        }

        R <- H - Omega%*%EG

        err <- colSums(R^2)/((n-m)*sum(X^2))
        ##rerr <- crossprod(X^2, R^2)/sum(X^2)^2 ##robust standard errors HC0
        se <- as.vector(sqrt(err)/sqrt(sigma2))
        beta <- as.vector(EG[1,])
        t <- beta/se
        pr <- 2*pnorm(-abs(t))
        list(loglik=loglik[1:i], beta=beta, se=se, t=t, p=p, R = R, loadings=Omega[,(q+2):m], factors=EG[(q+2):m,])
    }

    lrt.heft <- function(Y, Z, X, k, lambda=1, niter=50){
        loglik.full <- heft(Y, Z, X, k, lambda, niter)$loglik
        loglik.full <- loglik.full[length(loglik.full)]
        p <- ncol(Y)
        LR <- numeric(100)
        for(j in 1:100){
            loglik <- heft(Y[,-j], Z, X, k, lambda, niter)$loglik
            LR[j] <- -2*(loglik[length(loglik)] - loglik.full)
        }
        list(lr=LR, p=pchisq(LR, 1))
    }

    id <- covariates[,7] == 3

    ##Age
    Y <- log2(t(counts[,id])+2)
    Z <- covariates[id,-c(4, 7, 11)]
    X <- matrix(covariates[id, 4], ncol=1)

    ##smoking
    Y <- log2(t(counts[,id])+2)
    Z <- covariates[id,-c(6, 7, 11)]
    X <- matrix(covariates[id, 6], ncol=1)

    pnull <- matrix(0, ncol(Y), 100)
    for(i in 1:100) {
        Xp <- X[sample(1:nrow(X)),1, drop=FALSE]
        pnull[,i] <- heft(Y, Z, Xp, 1)$p
    }

    hist(as.vector(pnull), n=100)

    library(car)
    js <- sample(ncol(Y), 1000)
    lams <- as.numeric(1000)
    for(j in 1:1000)
        lams[j] <- powerTransform(Y[, js[j]])$roundlam

    YT <- bcPower(Y, lams)

    str(YT)

    plot(h$loglik, type='b', xlab="iteration", ylab="loglik", bty='n', lwd=2, pch=16, cex=1.6)
    grid(col=1)

    library(qqman)
    p <- h$p
    qq(p)
    names(p) <- colnames(Y)
    padj <- p.adjust(p, method="bonferroni")
    padj <- sort(padj)
    genes <- names(head(padj, n=25))
    lambda(p)

    plot(h$beta, h$se, col = as.numeric(p < 0.05/ncol(Y))+1, xlab="beta", ylab="SE(beta)")
    grid()
    sum((p < 0.05/ncol(Y)))

    library(org.Hs.eg.db)
    library(AnnotationDbi)
    AnnotationDbi::select(org.Hs.eg.db, keys=genes, keytype="ENSEMBL", column="SYMBOL")

    library(MASS)
    gene <- which(colnames(Y) == genes[1])
    cor(X, Y[, gene], method="spearman")
    plot(X, Y[, gene], ylab="Expression (log counts per million)", xlab="Age (Years)")
    abline(lm(Y[, gene]~X), col=2)
    abline(rlm(Y[, gene]~X), col=3)
    
}
