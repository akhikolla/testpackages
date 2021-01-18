nroMatch <- function(
    centroids,
    data) {

    # Check if input is a list.
    som <- list()
    if(!is.data.frame(centroids) && is.list(centroids)) {
	som <- centroids
        centroids <- centroids$centroids
    }
    if(length(centroids) < 1) {
        warning("Empty input.")
        return(NULL)
    }

    # Check variable names.
    vars <- colnames(centroids)
    if(length(vars) < 1) stop("No column names.")
    vars <- intersect(vars, colnames(data))
    if(length(vars) < 1) stop("Incompatible inputs.")
    if(length(vars) < ncol(centroids))
        warning("Incomplete coverage of variables.")

    # Convert inputs to numeric matrices.
    centroids <- nroRcppMatrix(centroids[,vars], trim=FALSE)
    data <- nroRcppMatrix(data[,vars], trim=FALSE)

    # Find best-matching units.
    res <- .Call("nro_match",
        as.matrix(centroids),
        as.matrix(data),
        PACKAGE="Numero")
    if(is.character(res)) stop(res)
    
    # Convert to data frame.
    res <- data.frame(res, stringsAsFactors=FALSE)

    # Check if training history is available.
    delta <- NA; sigma <- NA
    if(is.null(som$history) == FALSE)
        delta <- som$history[length(som$history)]
    if(is.null(som$layout) == FALSE) {
        sigma <- stats::quantile(som$layout$RESIDUAL,
	    c(0.3085, 0.6915), na.rm=TRUE)
        sigma <- (sigma[2] - sigma[1])
    }

    # Set mismatched labels to NA.
    bmus <- as.integer(res$DISTRICT)
    bmus[which(bmus == 0)] <- NA
    res$DISTRICT <- NULL

    # Standardize residuals against the the training error.
    res$RESIDUAL.z <- (res$RESIDUAL - delta)/(sigma + 1e-9)

    # Separate primary output from other information.
    names(bmus) <- rownames(data)
    rownames(res) <- rownames(data)
    attr(bmus, "quality") <- res
    attr(bmus, "variables") <- vars
    return(bmus)
}
