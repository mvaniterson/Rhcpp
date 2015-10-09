##convert seconds to minutes, hours and days
.sec2time <- function(x) {
    time <- if(x < 60)
        paste0(round(x, 0), "s.")
    else if(x < 60*60)
        paste0(x%/%60, "m:", .sec2time(x%%60))
    else if(x < 60*60*24)
        paste0(x%/%(60*60), "h:", .sec2time(x%%(60*60)))
    else if(x < 60*60*24*7)
        paste0(x%/%(60*60*24), "d:", .sec2time(x%%(60*60*24)))
    time
}

##draw samples from an inverse Gamma-distribution
.rigamma <- function(n, a, b) {
    return(1/rgamma(n, shape=a, rate=b))
}

##draw samples from a Dirichlet-distribution
.rdirichlet <- function(n, par) {
    k <- length(par)
    z <- array(0, dim =c(n, k))
    s <- array(0, dim =c(n, 1))
    for(i in 1:k){
        z[,i] <- rgamma(n, shape = par[i])
        s <- s + z[,i]
    }
    for(i in 1:k) {
        z[,i] <- z[,i]/s
    }
    return(z)
}

