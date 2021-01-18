.emve.Rcpp <- function(x, x_nonmiss, pu, n, p, theta0, G, d, x.miss.group.match, miss.group.unique, miss.group.counts,
        miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, nresample, nsubsize, sampling, cc, ck, EM_maxits){

    x_na <- x
    x_na[x_nonmiss == 0] <- NA
    subsample.id <- 0:(n-1)
    if(sampling == "cluster") subsample.id <- .sample.cluster(x_na) - 1
    
    res <- tryCatch( .Call("emve_Rcpp", x, x_nonmiss, pu, n, p, theta0, G, d, x.miss.group.match, miss.group.unique, miss.group.counts,
        miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, nresample, nsubsize, subsample.id, cc, ck, EM_maxits),
        "std::range_error" = function(e){
        conditionMessage( e ) } )

    if( is.character( res ) )
        stop( paste(res, "\nEstimated covariance matrix is not positive definite.\n") )

    ss0 <- res[1,1,]
    mu0 <- t(res[2,,])
    S0 <- res[-c(1,2),,]

    list(ss0=ss0, mu0=mu0, S0=S0)
}

###################################################################
## Compute the correction constants k and 50% chi-sq quantiles for 
## the computation of the scale of MVE for incomplete data (new version)
.scale.mve.init <- function(pp){
    cc <- qchisq(1/2,pp)
    k1 <- 1/pp
    k2 <- 0.5^(pp/2)
    k3 <- 1/(gamma(pp/2))
    k4 <- cc^(1+pp/2)
    k5 <- exp(-cc/2)
    kk <- k1 * k2 * k3 * k4 * k5
    ck <- cc*kk
    return( list(cc=cc, ck=ck) )
}


###################################################################
## Partial mahalanobis distance
## Rcpp version
.partial.mahalanobis.Rcpp <- function(x_mu_diff, Sigma, miss_group_unique, miss_group_counts){
    res <- tryCatch( .Call("fast_partial_mahalanobis", x_mu_diff, Sigma, miss_group_unique, miss_group_counts),
        "std::range_error" = function(e){
        conditionMessage( e ) } )
    if( is.character( res ) ) stop(res)
    return( c(res) )
}


###################################################################
## clustering-based sampling
.sample.cluster <- function(x, explvar = 0.9999, method.clust="ward.D", scalefn=Qn)
{
    n <- nrow(x)
    p <- ncol(x)
    x.loc <- apply(x, 2, median, na.rm=TRUE)
    x.scale <- apply(x, 2, function(x) scalefn(na.omit(x)))
    if(any(x.scale == 0)) stop("More than 50% equal values in one or more variables!")
    x.sc <- scale(x, x.loc, x.scale)
    R <- matrix(1, p, p)
    for (i in 1:(p-1)) for (j in (i+1):p) {
        R[j, i] <- R[i, j] <- (scalefn( na.omit(x.sc[, i] + x.sc[, j]) )^2 - 
                                 scalefn( na.omit(x.sc[, i] - x.sc[, j]) )^2)/4
    }
    if(any(is.na(x))) x.sc <- .impute.coord.med(x)
    R.eigen <- eigen(R)
    R.eigen.L <- R.eigen$values
    p0 <- rev(which(R.eigen.L > 0))[1]
    R.eigen.L[which(R.eigen.L < 0)] <- 0
    p1 <- (1:p0)[(cumsum(R.eigen.L)/sum(R.eigen.L) > explvar)][1]
    x.pc <- x.sc %*% R.eigen$vectors[, 1:p1]
    xpc.sc <- scale(x.pc, apply(x.pc, 2, median), apply(x.pc, 2, scalefn))
    hctree <- hclust(dist(xpc.sc), method=method.clust)
    cont <- TRUE
    b <- 1
    while(cont){
        b <- b + 1
        clusterB <- cutree(hctree, k=b)
        cont <- !all(table(clusterB) < n/2)
    }
    b <- b-1
    clusterB <- cutree(hctree, k=b)
    clusterB_sizes <- table(clusterB)
    clusterB1 <- which( clusterB_sizes == max(clusterB_sizes))
    clusterB1.id <- which(clusterB %in% clusterB1)

    return(clusterB1.id)  
}
