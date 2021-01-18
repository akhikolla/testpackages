nroPermute <- function(
    map,
    districts,
    data,
    n=1000,
    message=NULL) {

    # Convert data to numeric matrix.
    data <- nroRcppMatrix(data, trim=FALSE)
    topology <- nroRcppMatrix(map$topology, trim=FALSE)

    # Check training variables.
    trvars <- attr(districts, "variables")
    if(length(trvars) < 1) trvars <- colnames(map$centroids)
    if(length(trvars) < 1)
        warning("No information on training variables.")
    if(length(colnames(data)) < 1)
        warning("No column names in input data.")
    trvars <- intersect(colnames(data), trvars)

    # Check districts.
    districts <- nroRcppVector(districts,
        default=NULL, numeric=is.numeric(districts))

    # Check smoothness.
    smoothness <- attr(map$topology, "smoothness")
    smoothness <- nroRcppVector(smoothness[[1]], default=1)
    if(!is.finite(smoothness)) stop("Unusable map smoothness.")
    if(smoothness < 1) stop("Map smoothness less than one.")

    # Check parameters.
    n <- as.integer(nroRcppVector(n[[1]], default=1000))
    message <- nroRcppVector(message[[1]], default=-1)

    # Remove empty data columns.
    mu <- colMeans(data, na.rm=TRUE)
    cols <- which(0*mu == 0)
    if(length(cols) < length(mu)) {
        warning("Unusable columns excluded.")
        data <- data[,cols]
    }

    # Check input sizes.
    if(nrow(data) != length(districts))
        stop("Incompatible inputs.")

    # Check maximum number of cycles.
    if(!is.finite(n)) stop("Unusable number of permutations.")
    if(n < 10) stop("Too few permutations.")

    # Check message interval.
    if(!is.finite(message)) message <- -1.0

    # Collect training samples.
    trkeys <- intersect(rownames(data), names(districts))
    if(length(rownames(data)) < 1)
        warning("No row names in input data.")
    if(length(names(districts)) < 1)
        warning("No identity information in district data.")

    # Set flags for training variables.
    trflags <- rep("no", length.out=ncol(data))
    cols <- which(match(colnames(data), trvars) > 0)
    if(length(trkeys) > 0) trflags[cols] <- "yes"
    if(length(rownames(data)) < 1) trflags[cols] <- "unknown"
    if(length(names(districts)) < 1) trflags[cols] <- "unknown"
    if(ncol(data) == 1) trflags <- "unknown"

    # Set maximum number of permutations.
    evmask <- which(trflags != "yes")
    trmask <- which(trflags == "yes")
    numcycl <- rep(n, length.out=ncol(data))
    numcycl[trmask] <- min(1000, n) 

    # Estimate statistics.
    res <- .Call("nro_permute",
        as.matrix(topology),
        as.double(smoothness),
        as.integer(districts),
        as.matrix(data),
        as.integer(numcycl),
        as.double(message),
        PACKAGE="Numero") 
    if(is.character(res)) stop(res)

    # Convert results to data frame.
    output <- data.frame(res, stringsAsFactors=FALSE)
    rownames(output) <- colnames(data)
    output$TRAINING <- trflags

    # Estimate P-values.
    output$P.z <- stats::pnorm(output$Z, lower.tail=FALSE)
    output$P.z[trmask] <- NA
    output$P.freq[trmask] <- NA
    output$AMPLITUDE <- NA

    # Calculate base scores.
    z.tr <- output$Z[trmask]
    z.ev <- output$Z[evmask]
    trbase <- stats::quantile(z.tr, probs=0.95, na.rm=TRUE)
    evbase <- stats::quantile(z.ev, probs=0.95, na.rm=TRUE)
    if(!is.finite(trbase)) trbase <- evbase
    if(!is.finite(evbase)) evbase <- trbase
    if(!is.finite(trbase)) {
        warning("Permutations failed.")
        return(output)
    }

    # Attenuate high-scoring training variables.
    if((trbase > evbase) && (evbase > 1)) {
        mask <- which(z.tr > evbase)
	delta <- 0.2*sqrt(z.tr[mask] - evbase)
        z.tr[mask] <- (evbase + delta)
    }

    # Estimate color amplitudes.
    z <- c(z.tr, z.ev)
    rows <- c(trmask, evmask)
    zbase <- stats::quantile(z, probs=0.95, na.rm=TRUE)
    zbase <- max(1.1*zbase, 3)
    output$AMPLITUDE[rows] <- pmax(z/zbase, 0.04)
    attr(output, "zbase") <- zbase
    return(output)
}
