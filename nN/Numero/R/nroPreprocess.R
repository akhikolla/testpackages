nroPreprocess <- function(
    data,
    method="standard",
    clip=3.0,
    resolution=100,
    trim=FALSE) {

    # Convert input to numeric matrix.
    data <- nroRcppMatrix(data, trim=trim[[1]])
    binary <- attr(data, "binary")

    # Check if any rows or columns were excluded.
    if(length(attr(data, "excl.rows")) > 0)
        warning("Unusable rows excluded.")
    if(length(attr(data, "excl.columns")) > 0)
        warning("Unusable columns excluded.")

    # Check input size.
    if(nrow(data) < 1) {
        warning("No usable data.")
        return(NULL)
    }

    # Check method.
    method <- as.character(method[[1]])
 
    # Check resolution.
    resolution <- as.integer(resolution[[1]])
    if(resolution < 20) # see downsampling
        stop("Unusable resolution.")

    # Standardize location and scale.
    ds.in <- data
    ds.out <- NA*ds.in
    for(vn in colnames(ds.out))
        ds.out[,vn] <- nroPreprocess.std(ds.in[,vn], method, clip)

    # Downsample data model.
    model <- nroPreprocess.down(ds.in, ds.out, resolution, method)

    # Truncate extreme values.
    for(vn in colnames(ds.out))
        ds.out[,vn] <- nroPreprocess.clip(ds.out[,vn], method, clip)

    # If no preprocessing, binary variables remain binary.
    if(method == "") {
        binary <- intersect(binary, colnames(ds.out))
        attr(ds.out, "binary") <- binary
    }

    # Return results.
    attr(ds.out, "mapping") <- model
    return(ds.out)
}

#---------------------------------------------------------------------------

nroPreprocess.std <- function(x, method, clip=NA) {

    # Check variance.
    sigma <- stats::sd(x, na.rm=TRUE)
    if(!is.finite(sigma)) return(x)
    if(sigma <= .Machine$double.eps) return(x)

    # No standardization.
    if(length(method) < 1) return(x)
    if(nchar(method) < 1) return(x)

    # Rank-based standardization.
    if((method == "uniform") || (method == "tapered")) {
        z <- rank(x, na.last="keep")
        z <- (z - min(z, na.rm=TRUE))
	z <- (2*z/max(z, na.rm=TRUE) - 1)
	if(method == "tapered") z <- (z + 2*z^3)/3
        return(z)
    }

    # Default method left.
    if(method != "standard") stop("Unknown method.")

    # Protect against extreme outliers.
    t <- nroPreprocess.clip(x, method="standard", clip=5.0)
    t <- stats::na.omit(t)

    # Check if logarithm is useful.
    tmin <- min(t, na.rm=TRUE)
    if((tmin >= 0) && (sum(is.finite(t)) >= 10)) {
         t.log <- log(t + 1e-20)

         # Downsample for Shapiro test.
         mask <- which(0*t.log == 0)
         if(length(mask) > 5000)
	     mask <- sample(mask, size=5000)    

         # Test for normality.
         suppressWarnings(w <- stats::shapiro.test(t[mask]))
         suppressWarnings(w.log <- stats::shapiro.test(t.log[mask]))
	 if((w$p.value < 0.05) && (w$statistic < w.log$statistic)) {
             x <- log(x + 1e-20)
             t <- t.log
	 }	 
    }

    # Basic statistics.
    mu <- mean(t, na.rm=TRUE)
    sigma <- stats::sd(t, na.rm=TRUE)

    # Standardize scale and location.
    z <- (x - mu)/max(sigma, 1e-20)
    return(z)
}

#---------------------------------------------------------------------------

nroPreprocess.clip <- function(x, method, clip) {
    if((method != "standard") && (method != "")) return(x)
    if(!is.finite(clip)) return(x)
    med <- stats::median(x, na.rm=TRUE)
    sigma <- stats::sd(x, na.rm=TRUE)
    xmin <- (med - clip*sigma)
    xmax <- (med + clip*sigma)
    x[which(x < xmin)] <- xmin
    x[which(x > xmax)] <- xmax
    return(x)
}
    
#---------------------------------------------------------------------------

nroPreprocess.down <- function(x, y, resol, method) {
    if(method == "") return(NULL)

    # Nothing to do.
    results <- list()
    results$input <- x
    results$output <- y
    if(nrow(x) <= resol) return(results)

    # Prepare result matrices.
    results$input <- matrix(NA, nrow=resol, ncol=ncol(x))
    results$output <- matrix(NA, nrow=resol, ncol=ncol(x))
    colnames(results$input) <- colnames(x)
    colnames(results$output) <- colnames(x)

    # Reduce resolution.
    ranked <- (method == "uniform") || (method == "tapered")
    for(vn in colnames(x)) {
       rows <- is.finite(x[,vn]*y[,vn])
       u <- x[rows,vn]
       v <- y[rows,vn]

       # Remove duplicates.
       mask <- which(!duplicated(u))
       u <- u[mask]
       v <- v[mask]
       n <- length(u)
       if(n < 2) next

       # Sort by input value.
       sorted <- order(u)
       u <- u[sorted]
       v <- v[sorted]

       # Set sentinel points.
       if(!ranked) {
           q <- c(1, 10, (resol - 10), (resol - 1))/resol
           sigma.u <- stats::quantile(u, c(0.01, 0.1, 0.9, 0.99), na.rm=T)
           sigma.v <- stats::quantile(v, c(0.01, 0.1, 0.9, 0.99), na.rm=T)
           delta.u <- diff(sigma.u)
           delta.v <- diff(sigma.v)
           u <- c((u[1] - 3*delta.u[1]), u, (u[n] + 3*delta.u[3]))
           v <- c((v[1] - 3*delta.v[1]), v, (v[n] + 3*delta.v[3]))
       }

       # Select sampling points.
       n <- length(u)
       pivots <- seq(from=2, to=(n-1), length.out=(resol-2))
       pivots <- c(1, pivots, n)
       u.pivots <- stats::approx(x=(1:n), y=u, xout=pivots)$y

       # Interpolate output values.
       results$input[,vn] <- u.pivots
       results$output[,vn] <- stats::approx(x=u, y=v, xout=u.pivots)$y
    }
    return(results)
}
