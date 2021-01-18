## code from Ben Goodrich posted here
## https://groups.google.com/forum/#!msg/stan-users/3gDvAs_qwN8/Xpgi2rPlx68J
##

rgbeta <- function(num, shape) {
    if(shape == Inf)     rep(0, num)
    else if(shape > 0)  -1 + 2 * rbeta(num, shape, shape)
    else if(shape == 0) -1 + 2 * rbinom(num, 1, 0.5)
    else stop("shape must be non-negative")
}

## LKJ
rcorvine <- function(n, eta = 1, cholesky = FALSE, permute = !cholesky) {
    alpha <- eta + (n - 2) / 2
    L <- matrix(0, n, n)
    L[1,1] <- 1
    L[-1,1] <- partials <- rgbeta(n - 1, alpha)
    if(n == 2) {
        L[2,2] <- sqrt(1 - L[2,1]^2)
        if(cholesky) return(L)
        Sigma <- tcrossprod(L)
        if(permute) {
            ord <- sample(n)
            Sigma <- Sigma[ord,ord]
        }
        return(Sigma)
    }
    W <- log(1 - partials^2)
    for(i in 2:(n - 1)) {
        gap <- (i+1):n
        gap1 <- i:(n-1)
        alpha <- alpha - 0.5
        partials <- rgbeta(n - i, alpha)
        L[i,i] <- exp(0.5 * W[i-1])
        L[gap,i] <- partials * exp(0.5 * W[gap1])
        W[gap1] <- W[gap1] + log(1 - partials^2)
    }
    L[n,n] <- exp(0.5 * W[n-1])
    if(cholesky) return(L)
    Sigma <- tcrossprod(L)
    if(permute) {
        ord <- sample(n)
        Sigma <- Sigma[ord,ord]
    }
    return(Sigma)
}

