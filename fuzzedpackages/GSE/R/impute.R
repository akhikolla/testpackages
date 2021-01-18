.impute.coord.med <- function(x){
    x.med <- apply(x, 2, median, na.rm=TRUE)
    u <- is.na(x)
    x[is.na(x)] <- 0
    x <- x + sweep(u, 2, x.med, "*")
    x
}

ImpS <- function(x, alpha=0.95, method=c("bisquare","rocke"), init=c("emve","emve_c"), ...){
    xcall <- match.call()
    
    method <- match.arg(method)
    init <- match.arg(init)
    
    impute <- .impute.coord.med

    x.filt <- gy.filt(x, alpha=c(alpha, 0))
    x.imp <- impute(x.filt)
    res <- GSE(x.imp, method=method, init=init, ...)

    res <- new("TSGS",
        call = xcall,
        S = res@S,
        mu = res@mu,
        xf = x.filt,
        sc = res@sc,
        mu0 = res@mu0,
        S0 = res@S0, 
        iter = res@iter,
        eps = res@eps,
        estimator = "2SGS", 
        x = x,
        ximp = res@ximp,
        weights = res@weights,
        weightsp = res@weightsp,
        pmd = res@pmd,
        pmd.adj = res@pmd.adj,
        p = res@p,
        pu = res@pu)
        
    return(res)
}
