numero.evaluate <- function(
    model,
    data,
    ranked=TRUE,
    n=1000) {

    # Continue analyses.
    output <- list(stamp=date())
    cat("\n*** numero.evaluate ***\n", output$stamp, "\n", sep="")

    # Check that resources are available.
    if(is.null(model$map)) stop("Self-organizing map not available.")
    if(is.null(model$layout)) stop("Map layout not available.")
    layout <- model$layout

    # Check if input is a vector.
    if(is.vector(data)) stop("Data must be a matrix or a data frame.")

    # Check that data and layout are compatible.
    cat("\nDataset:\n")
    pos <- match(rownames(data), rownames(layout))
    rows <- which(pos > 0)
    if(length(rows) < 1) {
        warning("Incompatible data and layout.")
        return(NULL)
    }

    # Harmonize data and layout.
    nprev <- nrow(data)
    data <- data[rows,,drop=FALSE]
    layout <- layout[pos[rows],,drop=FALSE]
    cat(nrow(data), " / ", nprev, " rows included\n", sep="")
    cat(ncol(data), " columns included\n", sep="")

    # Add identities to district assignments.
    bmc <- layout[,"BMC"]
    names(bmc) <- rownames(layout)

    # Apply rank transform to protect from extreme values.
    if(ranked) data <- nroPreprocess(data=data, method="uniform")

    # Calculate component planes.
    comps <- nroAggregate(topology=model$map, districts=bmc, data=data)

    # Estimate statistics in chunks.
    cat("\nStatistics:\n")
    suppressWarnings(
        stats <- nroPermute(map=model$map, districts=bmc,
                            data=data, n=n, message=10))
    cat(nrow(stats), " usable variables\n", sep="")
    cat(sum(stats$N.cycles), " permutations\n", sep="")

    # Make sure all variables are included.
    missed <- setdiff(colnames(data), rownames(stats))
    if(length(missed) > 0) {
        x <- list()
        for(v in colnames(stats))
          x[[v]] <- rep(NA, length(missed))
        x <- data.frame(x, stringsAsFactors=FALSE)
	rownames(x) <- missed
	stats <- rbind(stats, x)
	stats <- stats[colnames(data),,drop=FALSE]
    }

    # Revert transform.
    if(ranked) comps <- nroPostprocess(data=comps,
        mapping=attr(data,"mapping"), reverse=TRUE)

    # Determine district ranges.
    colrs <- nroColorize(comps)

    # Collect results.
    output$map <- model$map
    output$layout <- layout
    output$planes <- comps
    output$ranges <- attr(colrs, "ranges")
    output$palette <- "rhodo"
    output$statistics <- stats
    output$data <- data
    return(output)
}

#----------------------------------------------------------------------------

numero.evaluate.transf <- function(x, method) {
    output <- NULL
    mu <- mean(x, trim=0.1, na.rm=TRUE)
    sigma <- mean(abs(x - mu), trim=0.1, na.rm=TRUE)
    if(!is.finite(sigma)) return(NULL)
    if(method[[1]] == "probit") {
        output <- (x - mu)/sigma
	posit <- which(output > 0)
	negat <- which(output < 0)
	output[posit] <- log(1 + output[posit])
	output[negat] <- -log(1 - output[negat])
    }
    if(method[[1]] == "logarithm")
        output <- log(1 + x/mu)
    attr(output, "method") <- method
    attr(output, "param") <- c(mu, sigma)
    return(output)
}

#----------------------------------------------------------------------------

numero.evaluate.itransf <- function(x, method, param) {
    mu <- param[1]
    sigma <- param[2]
    alpha <- param[3]
    if(method[[1]] == "probit") {
        posit <- which(x > 0)
	negat <- which(x < 0)
	x[posit] <- (exp(x[posit]) - 1)
	x[negat] <- (1 - exp(-x[negat]))
	output <- (x*sigma + mu)
    }
    if(method[[1]] == "logarithm")
        output <- (exp(x) - 1)*mu
    return(output)
}
