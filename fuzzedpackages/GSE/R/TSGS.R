
TSGS <- function(x, filter=c("UBF-DDC","UBF","DDC","UF"),
            partial.impute=FALSE, tol=1e-4, maxiter=150,
            method=c("bisquare","rocke"),
            init=c("emve","qc","huber","imputed","emve_c"), mu0, S0){

    xcall <- match.call()

    filter <- match.arg(filter)
    method <- match.arg(method)
    init <- match.arg(init)

    ## check dat
    if(is.data.frame(x) | is.matrix(x))
        x <- data.matrix(x)
    else stop("Data matrix must be of class matrix or data.frame.")
    if(any(is.na(x))) warning("Data matrix contains missing values.")
    ## June 12, 2012
    ## Only allow up to p=50
    n <- nrow(x)
    p <- ncol(x)
    if( p >200 | p < 2 ) stop("Column dimension of 'x' must be in between 2 and 200.")

    ## 1) filter step
    if(filter == "UF"){
      xf <- gy.filt(x, alpha = c(0.95, 0))
    }
    else if(filter == "UBF"){
      xf <- gy.filt(x)
    }
    else if(filter == "DDC"){
      tmp <- capture.output({res.DDC <- cellWise::DDC(x)})
  		xf <- x
  		xf[res.DDC$indcells] <- NA
    }
    else if(filter == "UBF-DDC"){
      # UBF
      xf.ubf <- gy.filt(x)
      v.ubf <- 1*is.na(xf.ubf)
      # DDC
      tmp <- capture.output({res.DDC <- cellWise::DDC(x)})
      xf.ddc <- x
      xf.ddc[res.DDC$indcells] <- NA
      v.ddc <- 1*is.na(xf.ddc)
      # combine the two filter results
      xf <- x
      xf[v.ubf == 1 & v.ddc == 1] <- NA
    }

    ## 2) partial imputation step as suggested in the rejoinder of Agostinelli et al (2015)
    xf_pi <- xf
    if( partial.impute ){
        ximp <- .impute.coord.med(xf_pi)
        aid <- which(rowSums(!is.na(xf_pi)) == p)
        uid <- which(rowSums(!is.na(xf_pi)) < p)
        n0 <- n/2 + (p+1)
        if( n0 > length(aid) ){
            fid <- sample( uid, n0 - length(aid))
            xf_pi[fid,] <- ximp[fid,]
        }
    }

    ## 3) estimation step
    res <- GSE(xf_pi, tol=tol, maxiter=maxiter, method=method, init=init, mu0, S0)

    res <- new("TSGS",
        call = xcall,
        S = res@S,
        mu = res@mu,
        xf = xf,
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
    res
}
