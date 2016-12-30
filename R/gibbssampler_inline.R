

gibbssampler <- function(y, niter=5000,
                         alpha=1.28, beta=0.36*var(y),
                         lam = c(0, 1, -1), tau =  c(1000, rep(100, 2)),
                         g = c(90,5,5)) {

    k <- 3
    n <- length(y)
    mu <- lam
    sig <- rep(sd(y)/k, k)
    p <- g/sum(g)
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

        gibbsmu[i,] <- rnorm(k, mean = (lam*tau + sj)/(nj + tau), sd = sqrt(mixparam$sig^2/(tau+nj)))

        mixparam$mu <- gibbsmu[i,]

        for(j in 1:k) {
            sj2[j] <- sum(as.numeric(z == j) * (y - mixparam$mu[j])^2)
        }

        gibbssig[i, ] <- sqrt(.rigamma(k, alpha + 0.5*(nj + 1), beta + 0.5*tau *(mixparam$mu -lam)^2 + 0.5*sj2))
        mixparam$sig <- gibbssig[i,]

        gibbsp[i,] <- .rdirichlet(1, par = nj + g)
        mixparam$p <- gibbsp[i,]

    }
    gs <- data.frame(p = gibbsp, mu = gibbsmu, sigma = gibbssig)
}

library(inline)




##draw samples from a Dirichlet-distribution
.rdirichlet <- function(n, g) {
    NumericMatrix  = z n, 3
    NumericVector = s n
    
    for(int j = 0; j<3; j++){
        NumericVector zj = rgamma(n, shape = g[j]);
        for(int i = 0; i<n; i++)
            s[i] = s[i] + z[i,j];
    }
    
    for(int i = 0; i<n; i++)
        for(int j = 0; j<3; j++)
            z[i,i] = z[i,j]/s;
    return z;
}




src <- '
Rcpp::NumericVector vy(y);
int niter = 5000;

alpha = 1.28;
beta = 0.36*var(y);
lam[0] = 0
lam[1] = 1
lam[2] = -1;
tau[0] = 1000
tau[1] = tau[2] = 100;
g[0] = 90
g[1] = g[2] = 5;

int k = 3;
int n = vy.size();
mu = lam;
sig[0] = sig[1] = sig[2] = sd(y)/3;
p[0] = g[0]/sum(g);
p[1] = p[2] = g[1]/sum(g);
nj = sj = sj2 = 0;
NumericVector z

gibbsmu = gibbssig = gibbsp = matrix(0, nrow= niter, ncol = k);

for (int i = 0; i < niter; i++) {

    for(int j = 0; j < n; j++) {
        NumericVector prob = dnorm(y[j], mu, sig);
        z[j] = sample(k, 1, prob=prob);
    }

    for (int j = 0; j < k; j++) {
        nj[j] <- sum(z == j+1);
        sj[j] <- sum(as.numeric(z == j+1)*y);
    }

    mu = gibbsmu[i,] = rnorm(k, (lam*tau + sj)/(nj + tau), sqrt(sig^2/(tau+nj)));

    for (int j = 0; j < k; j++)
        sj2[j] <- sum(as.numeric(z == j) * (y - mu[j])^2);

    sig = gibbssig[i, ] = sqrt(1/rgamma(k, alpha + 0.5*(nj + 1), beta + 0.5*tau *(mu -lam)^2 + 0.5*sj2))

    p = gibbsp[i,] = .rdirichlet(1, par = nj + g);
}
gs <- data.frame(p = gibbsp, mu = gibbsmu, sigma = gibbssig)

return gs;
'

gibbssampler <- cxxfunction(signature(y="vector"), src, plugin = "Rcpp")
