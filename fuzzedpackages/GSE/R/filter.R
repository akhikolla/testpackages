.gy.filt.uni <- function(v, alpha, df){
    n <- length(v)
    id <- (1:n)[!is.na(v)]
    v.out <- rep(NA, n)
    
    v <- na.omit(v)
    n <- length(v)
    v.order <- order(v)	
    v <- sort(v)
    i0 <- which(v < qchisq( alpha, df ))
    n0 <- 0
    if( length(i0) > 0){
        i0 <- rev(i0)[1]
        dn <- max( pmax( pchisq( v[i0:n], df) - (i0:n - 1)/n, 0)) 
        n0 <- round(dn*n)
    } 
    v <- v[ order(v.order) ] 
    v.na <- v
    if(n0 > 0) v.na[ v.order[ (n - n0 + 1):n] ] <- NA
    
    v.out[id] <- v.na
    
    return( v.out )
}

.gy.filt.uni.it <- function(v, alpha, df, miter){
    converge <- 0
    iter <- 0
    n <- length(v)
    n_flag <- sum(is.na(v))
    while( converge == 0 & iter < miter ){
        iter <- iter + 1
        v <- .gy.filt.uni( v, alpha, df )
        n_flag_new <- sum(is.na(v))
        if( n_flag_new - n_flag <= 0 ) converge <- 1
    }
    return( v)
}

gy.filt <- function(x, alpha=c(0.95,0.85), bivarQt=0.99, bivarCellPr=0.1, miter=5){

    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    xs <- scale(x, apply(x, 2, median, na.rm=TRUE), apply(x, 2, mad, na.rm=TRUE))
    u <- 1-is.na(x)
    if( alpha[1] > 0) u <- 1-is.na(apply(xs^2, 2, .gy.filt.uni.it, alpha=alpha[1], df=1, miter=miter))
    m <- matrix(0, n, p)
    mn <- matrix(0, n, p)

    if( alpha[2] > 0){
        for(j in 1:(p-1))
            for(k in (j+1):p){
                S <- matrix(1, 2, 2)
                S[1,2] <- S[2,1] <- (mad( na.omit(xs[,j] + xs[,k]) )^2 - 
                                            mad( na.omit(xs[,j] - xs[,k]) )^2)/4

                pmd <- mahalanobis(xs[,c(j,k)], rep(0,2), S)
                pmd[ which(u[,j]*u[,k] == 0) ] <- -1
                pmd <- pmd * qchisq(0.5, df=2) / median(pmd[which(pmd > 0)])
                mn[,c(j,k)] <- mn[,c(j,k)] + 1*(pmd > 0)
                pmd[which(pmd > 0)] <- .gy.filt.uni.it(pmd[which(pmd > 0)], alpha=alpha[2], df=2, miter=miter)
                m[,c(j,k)] <- m[,c(j,k)] + 1*is.na(pmd)
            }
        u[m > qbinom(bivarQt, mn, prob=bivarCellPr) & mn > 0] <- 0
    }
    xna <- x
    xna[u == 0] <- NA
    return( xna )
}