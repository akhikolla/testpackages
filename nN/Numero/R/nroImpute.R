nroImpute <- function(
    data,
    subsample=500,
    standard=TRUE,
    message=NULL) {

    # Convert input to numeric matrix.
    data <- nroRcppMatrix(data, trim=FALSE)
    binary <- attr(data, "binary")
    binary <- nroRcppVector(binary, default=NULL,
        numeric=is.numeric(binary))

    # Check input size.
    if(nrow(data)*ncol(data) < 1) {
        warning("Empty input.")
        return(NULL)
    }
    if(nrow(data) < 2) return(data)
    if(ncol(data) < 2) return(data)

    # Ensure inputs are safe for C++.
    subsample <- nroRcppVector(subsample[[1]], default=nrow(data))
    standard <- (nroRcppVector(standard[[1]], default=1) == 1)
    message <- nroRcppVector(message[[1]], default=-1)

    # Check subsample size.
    if(!is.finite(subsample)) stop("Unusable subsample size.")
    if(subsample > nrow(data)) subsample <- nrow(data)
    if(subsample < 10) stop("Too small subsample.")

    # Check message interval.
    if(!is.finite(message)) stop("Unusable message interval.")

    # Copy names.
    rnames <- rownames(data)
    cnames <- colnames(data)
    
    # Detect numeric variables.
    numerics <- c()
    for(j in 1:ncol(data)) {
	flags <- is.finite(data[,j])
        if(sum(flags) > 0) numerics <- c(numerics, j)
    }
    if(length(numerics) < 2) {
        warning("Less than two numeric columns.")
        return(data)
    }
    if(length(numerics) < ncol(data))
        warning("Non-numeric columns.")

    # Standardize data.
    sigma <- rep(1, ncol(data))
    if(standard) {
        for(j in numerics) {
	    x <- data[,j]
            s <- stats::sd(x, na.rm=TRUE)
            if(!is.finite(s)) next
            if(s <= 0.0) next
	    sigma[j] <- s
            data[,j] <- x/s
	}
    }

    # Impute missing values.
    res <- .Call("nro_impute",
        as.matrix(data[,numerics,drop=FALSE]),
        as.integer(subsample),
        as.double(message),
        PACKAGE="Numero")
    if(is.character(res)) stop(res)
    data[,numerics] <- res

    # Restore original scale.
    for(j in numerics)
        data[,j] <- (sigma[j])*(data[,j])

    # Restore names and attributes.
    rownames(data) <- rnames
    colnames(data) <- cnames
    attr(data, "binary") <- intersect(binary, cnames)
    return(data)
}
