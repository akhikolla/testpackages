nroLabel <- function(
    topology,
    values,
    gap=2.3) {

    # Check if input is a list.
    if(!is.data.frame(topology) && is.list(topology))
        topology <- topology$topology

    # Convert inputs to numeric matrices.
    topology <- nroRcppMatrix(topology, trim=FALSE)
    values <- nroRcppMatrix(values, trim=FALSE)
    if(nrow(values)*ncol(values) < 1) {
        warning("Empty input.")
	return(NULL)
    }

    # Check topology and values.
    if(nrow(topology) < 2) stop("Unusable topology.")
    if(nrow(topology) != nrow(values)) stop("Incompatible inputs.")
	
    # Check gap.
    gap <- nroRcppVector(gap[[1]], default=2.3)
    if(!is.finite(gap)) stop("Unusable gap.")
    if(gap < 1.0) stop("Gap is less than one.")

    # Set flags for binary data.
    binflags <- match(colnames(values), attr(values, "binary"))
    binflags <- is.finite(binflags)

    # Determine label positions.
    res <- .Call("nro_label",
        as.matrix(topology),
        as.matrix(values),
	as.integer(binflags),
        as.numeric(gap),
        PACKAGE = "Numero" )
    if(is.character(res)) stop(res)

    # Convert to matrices.
    labels <- matrix("", nrow=nrow(values), ncol=ncol(values))
    visible <- matrix(NA, nrow=nrow(values), ncol=ncol(values))
    for(j in 1:ncol(values)) {
        labels[,j] <- res$labels[[j]]
        visible[,j] <- res$visible[[j]]
    }

    # Set row and column names.
    rownames(labels) <- rownames(values)
    rownames(visible) <- rownames(values)
    colnames(labels) <- colnames(values)
    colnames(visible) <- colnames(values)

    # Return results.
    attr(labels, "visible") <- visible
    return(labels)
}
