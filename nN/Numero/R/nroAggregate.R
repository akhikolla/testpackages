nroAggregate <- function(
    topology,
    districts,
    data=NULL) {

    # Check if input is a list.
    if(!is.data.frame(topology) && is.list(topology))
        topology <- topology$topology

    # Check smoothness.
    smoothness <- attr(topology, "smoothness")
    smoothness <- nroRcppVector(smoothness[[1]], default=NA)
    if(!is.finite(smoothness)) stop("Unusable map smoothness.")
    if(smoothness < 0) stop("Negative map smoothness.")

    # Ensure topology is a numeric matrix.
    topology <- nroRcppMatrix(topology, trim=FALSE)

    # Check if data are available.
    if((length(districts) < 1) || (ncol(topology) < 1)) {
        warning("Empty input.")
	return(NULL)
    }

    # Estimate sample histogram.
    if(is.null(data)) {
        res <- .Call("nro_diffuse",
            as.matrix(topology),
            as.double(smoothness),
	    as.integer(districts),
            matrix(nrow=0, ncol=0),
            PACKAGE="Numero")
        if(is.character(res)) stop(res)
	return(as.numeric(res$histograms))
    }

    # Convert data to numeric matrix.
    data <- nroRcppMatrix(data, trim=FALSE)
    binary <- attr(data, "binary")

    # Flag non-empty columns.
    mu <- colMeans(data, na.rm=TRUE)
    empty <- which(!is.finite(mu))
    if(length(empty) == ncol(data)) {
       warning("No usable data.")
       return(NULL)
    }

    # Replace empty columns with zeros.
    data[,empty] <- 0

    # Check compatibility.
    if(nrow(data) != length(districts))
        stop("Incompatible inputs.")

    # Estimate component planes.
    res <- .Call("nro_diffuse",
                 as.matrix(topology),  
                 as.double(smoothness),
                 as.integer(districts),
                 as.matrix(data),
                 PACKAGE="Numero");
    if(is.character(res)) stop(res)

    # Transpose to column-major format.
    planes <- t(res$planes)
    hgrams <- t(res$histograms)

    # Set row and column names.
    colnames(planes) <- colnames(data)
    colnames(hgrams) <- colnames(data)
    rownames(planes) <- 1:nrow(planes)
    rownames(hgrams) <- 1:nrow(hgrams)

    # Clear empty variables.
    planes[,empty] <- NA
    hgrams[,empty] <- 0

    # Return results.
    attr(planes, "histogram") <- hgrams
    attr(planes, "binary") <- intersect(binary, colnames(data))
    return(planes)
}
