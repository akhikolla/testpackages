nroKohonen <- function(
    seeds,
    radius=3,
    smoothness=1.0) {

    # Check if input is a list.
    if(!is.data.frame(seeds) && is.list(seeds))
        seeds <- seeds$centroids

    # Check input data.
    seeds <- nroRcppMatrix(seeds, trim=TRUE)
    if((nrow(seeds) < 3) || (ncol(seeds) < 3)) {
        warning("Not enough usable rows or columns.")
        return(NULL)
    }
	 
    # Check if any rows or columns were excluded.
    if(length(attr(seeds, "excl.rows")) > 0)
        warning("Unusable rows excluded.")
    if(length(attr(seeds, "excl.columns")) > 0)
        warning("Unusable columns excluded.")

    # Check radius.
    radius <- nroRcppVector(radius[[1]], default=3)
    if(!is.finite(radius)) stop("Unusable radius.")
    if(radius < 2) stop("Radius is less than two.")

    # Check smoothness.
    smoothness <- nroRcppVector(smoothness[[1]], default=NA)
    if(!is.finite(smoothness)) stop("Unusable smoothness.")
    if(smoothness < 1) stop("Smoothness less than one.")
    if(smoothness > 0.49*radius) stop("Smoothness too high.")

    # Set up a self-organizing map.
    res <- .Call("nro_kohonen",
        as.matrix(seeds),
        as.integer(radius),
        as.double(smoothness),
        PACKAGE="Numero");
    if(is.character(res)) stop(res)

    # Convert to data frame to make it easier to add columns later.
    res$topology <- data.frame(res$topology, stringsAsFactors=FALSE)

    # Set column names.
    colnames(res$centroids) <- colnames(seeds)
    colnames(res$topology) <- c("X", "Y", "RADIUS1", "RADIUS2",
                                "ANGLE1", "ANGLE2")
  
    # Set row names.
    rownames(res$centroids) <- (1:nrow(res$centroids))
    rownames(res$topology) <- (1:nrow(res$topology))

    # Save topology parameters.
    t <- table(res$topology$RADIUS1)
    attr(res$topology, "radius") <- radius <- (length(t) - 1)
    attr(res$topology, "smoothness") <- smoothness
    return(res)
}
