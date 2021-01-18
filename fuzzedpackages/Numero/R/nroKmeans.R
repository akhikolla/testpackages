nroKmeans <- function(
    data,
    k=3,
    subsample=NULL,
    balance=0,
    message=NULL) {

    # Convert data to numeric matrix.
    data <- nroRcppMatrix(data, trim=TRUE)
    if((nrow(data) < 10) || (ncol(data) < 3)) {
        warning("Not enough usable rows or columns.")
	return(NULL)
    }
	 
    # Check if any rows or columns were excluded.
    if(length(attr(data, "excl.rows")) > 0)
        warning("Unusable rows excluded.")
    if(length(attr(data, "excl.columns")) > 0)
        warning("Unusable columns excluded.")

    # Check parameters.
    k <- as.integer(nroRcppVector(k[[1]]), default=3)
    balance <- nroRcppVector(balance[[1]], default=0)
    message <- nroRcppVector(message[[1]], default=-1)
    subsample <- as.integer(nroRcppVector(subsample[[1]], default=NA))

    # Automatic subsample.
    if(!is.finite(subsample)) {
        subsample <- 10*sqrt(nrow(data))*sqrt(k)
	subsample <- min(subsample, 0.95*nrow(data), na.rm=TRUE)
        if(subsample/nrow(data) < balance)
            subsample <- balance*nrow(data)
	subsample <- round(subsample)
    }

    # Check balance parameter.
    if(!is.finite(balance)) stop("Unusable balance parameter.")
    if(balance < 0.0) stop("Unusable balance parameter.")
    if(balance > 1.0) stop("Unusable balance parameter.")

    # Check subsample size.
    if(!is.finite(subsample)) stop("Unusable subsample.")
    if(subsample > nrow(data)) subsample <- nrow(data)
    if(subsample < (k + 10)) stop("Too small subsample.")
    
    # Check message interval.
    if(!is.finite(message)) message <- -1.0
    
    # Estimate centroids.
    res <- .Call("nro_train",
        as.matrix(k),
        as.double(0.0),
        as.matrix(NA),
        as.matrix(data),
        as.integer(subsample),
        as.double(balance),
	as.double(message),
        PACKAGE="Numero")
    if(is.character(res)) stop(res)
    
    # Recode missing unit labels.
    res$layout[which(res$layout <= 0)] <- NA
    res$layout <- data.frame(BMC=res$layout, RESIDUAL=res$residuals)
    rownames(res$layout) <- rownames(data)
    res$residuals <- NULL

    # Copy column names.
    colnames(res$centroids) <- colnames(data)
    res$subsample <- subsample
    return(res)
}
