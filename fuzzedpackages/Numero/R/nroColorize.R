nroColorize <- function(
    values,
    ranges=NULL,
    amplitudes=1,
    palette="rhodo") {

    # Convert vector input to matrix.
    values <- nroRcppMatrix(values, trim=FALSE)

    # Check that ranges are usable.
    if(length(ranges) > 0) {
        if(is.vector(ranges) && (length(ranges) == 2))
	    ranges <- t(as.matrix(ranges))
        if(nrow(ranges) != ncol(values)) {
            warning("Unusable ranges, reverting to defaults.")
            ranges <- NULL
        }
    }

    # Default value ranges.
    if(length(ranges) < 1) {
        mu <- apply(values, 2, stats::median, na.rm=TRUE)
        delta <- apply(values, 2, function(x) {
            m <- stats::median(x, na.rm=TRUE)
	    d <- stats::quantile(abs(x - m), 0.99, na.rm=TRUE)
            return(d)
        })
        ranges <- data.frame(MIN=(mu - delta), MAX=(mu + delta),
                             stringsAsFactors=FALSE)
        rownames(ranges) <- colnames(values)
    }

    # Check palette.
    palette <- as.character(palette[[1]])

    # Check if amplitudes is a data frame or a matrix.
    if(is.data.frame(amplitudes) || is.matrix(amplitudes))
         amplitudes <- amplitudes[,"AMPLITUDE"]
    amplitudes <- nroRcppVector(amplitudes, default=1)

    # Expand scalar amplitudes.
    if(length(amplitudes) < ncol(values)) {	
        if(length(amplitudes) != 1)
	     warning("Amplitudes replicated to compatible number.")
        amplitudes <- rep(amplitudes, length.out=ncol(values))
    }

    # Center by the middle of the range.
    mu <- 0.5*(ranges[,1] + ranges[,2])
    z <- sweep(values, 2, mu, `-`)

    # Standardize by the width of the range.
    delta <- 0.5*(ranges[,2] - ranges[,1] + 1e-9)
    z <- sweep(z, 2, delta, `/`)

    # Amplification factors.
    z <- sweep(z, 2, amplitudes[1:ncol(z)], `*`)
    z <- 0.5*(z + 1.0) # centered at 0.5
    
    # Set colors.
    res <- .Call("nro_colorize",
        as.matrix(z),
        as.character(palette),
        PACKAGE="Numero")
    if(is.character(res)) stop(res)

    # Convert to data frame.
    if(ncol(values) < 2) {
        res$colors <- data.frame(X=as.vector(res$colors[[1]]),
	    stringsAsFactors=FALSE)
        res$contrast <- data.frame(X=as.vector(res$contrast[[1]]),
	    stringsAsFactors=FALSE)
    }
    else {
        res$colors <- data.frame(res$colors, stringsAsFactors=FALSE)
        res$contrast <- data.frame(res$contrast, stringsAsFactors=FALSE)
    }

    # Set row and column names.
    rownames(res$colors) <- rownames(values)
    rownames(res$contrast) <- rownames(values)
    colnames(res$colors) <- colnames(values)
    colnames(res$contrast) <- colnames(values)

    # Return results.
    attr(res$colors, "contrast") <- res$contrast
    attr(res$colors, "ranges") <- ranges
    return(res$colors)
}
