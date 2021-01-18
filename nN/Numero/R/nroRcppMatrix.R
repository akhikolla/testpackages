nroRcppMatrix <- function(
    data,
    trim) {

    # Check if anything to do.
    if(length(data) < 1) return(matrix(nrow=0, ncol=0))
    binary <- attr(data, "binary")

    # Convert to matrix.
    if(is.data.frame(data))
        data <- as.matrix(data)
    if(is.atomic(data) && !is.matrix(data)) {
        data <- as.matrix(data)
 	colnames(data) <- "data"
    }
    if(is.list(data)) {
        warning("Unusable input.")
	return(NULL)
    }

    # Check if already fully numeric.
    formatted <- FALSE
    if(is.matrix(data)) {
        if(is.numeric(data)) formatted <- TRUE
        if(is.integer(data)) formatted <- TRUE
        if(is.logical(data)) formatted <- TRUE
    }

    # Row and column names.
    rnames <- rownames(data)
    cnames <- colnames(data)
    if(length(rnames) != nrow(data)) {
        rnames <- (1:nrow(data))
	rownames(data) <- rnames
    }
    if(length(cnames) != ncol(data)) {
        cnames <- (1:ncol(data))
        colnames(data) <- cnames
    }

    # Convert to numeric values.
    if(!formatted) {
        try(
	suppressWarnings( # data frame to matrix
        data <- apply(data, 2, function(x) {
            x <- as.numeric(x)
            x[which(!is.finite(x))] <- NA
            return(x)
        })), silent=TRUE)
        if(!is.numeric(data)) {
            warning("Unusable input.")
	    return(NULL)
        }

        # Restore names.
        rownames(data) <- rnames
        colnames(data) <- cnames
    }

    # Detect binary variables.
    if(is.null(binary)) {
        for(vn in colnames(data)) {
            x <- data[,vn]
            n <- sum(is.finite(x), na.rm=TRUE)
            n0 <- sum((x == 0), na.rm=TRUE)
            n1 <- sum((x == 1), na.rm=TRUE)
            if(n != (n0 + n1)) next
            binary <- c(binary, vn)
        }
    }

    # Remove unusable data points.
    if(trim) {
        nprev <- 0
        while(nrow(data)*ncol(data) != nprev) {
            nprev <- nrow(data)*ncol(data)
	    if(nprev < 1) break

            # Remove rows with no usable values.
	    if(ncol(data) > 1) {
                rows <- which(is.finite(rowMeans(data, na.rm=TRUE)))
		data <- data[rows,,drop=FALSE]
            }
	    else {
	        rows <- which(is.finite(data[,1]))
		data <- data[rows,,drop=FALSE]
	    }

            # Remove columns with no usable values.
	    if(nrow(data) > 1) {
                cols <- which(is.finite(colMeans(data, na.rm=TRUE)))
		data <- data[,cols,drop=FALSE]
	    }
	    else {
                cols <- which(is.finite(data[1,]))
		data <- data[,cols,drop=FALSE]
            }
       }
    }

    # Return results.
    attr(data, "binary") <- intersect(binary, colnames(data))
    attr(data, "excl.rows") <- setdiff(rnames, rownames(data))
    attr(data, "excl.columns") <- setdiff(cnames, colnames(data))
    return(data)
}
