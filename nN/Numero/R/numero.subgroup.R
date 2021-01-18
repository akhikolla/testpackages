numero.subgroup <- function(
    results,
    variables,
    topology=NULL,
    reference=NULL,
    gain=1.0,
    detach=FALSE,
    capacity=9) {

    # Start processing.
    cat("\n*** numero.subgroup ***\n", date(), "\n", sep="")

    # Default inputs.
    if(is.null(variables)) variables <- colnames(results$planes)
    if(is.null(topology)) topology <- results$som$topology
    if(is.null(reference)) reference <- results
    detach <- as.character(detach[[1]])

    # Check capacity.
    cat("\nResources:\n")
    if(capacity < 2) {
        cat("capacity less than two\n")
        return(0)
    }

    # Select variables.
    comps <- results$planes
    stats <- results$statistics
    variables <- intersect(variables, colnames(comps))
    cat(length(variables), " columns included\n", sep="")
    if(length(variables) < 2) {
        cat("less than two usable variables\n")
        return(0)
    }

    # Check if too many variables.
    comps <- comps[,variables]
    stats <- stats[variables,]
    if(nrow(stats) > capacity) {
        cat("capacity exceeded, showing", capacity, "plots.\n")
        comps <- comps[,1:capacity]
	stats <- stats[1:capacity,]
    }

    # Check if reference is usable.
    rvars <- rownames(reference$statistics)
    if(sum(is.na(match(variables, rvars))) > 0) {
	cat("incompatible reference\n")
        return(0)
    }

    # Check if gain is usable.
    gain <- as.double(gain[[1]])
    if(!is.finite(gain)) {
        gain <- 1
        cat("unusable gain, reverted to 1\n")
    }
    if(gain <= 0.0) {
        gain <- 1
        cat("non-positive gain, reverted to 1\n")
    }

    # Get coloring parameters.
    amplitudes <- reference$statistics[variables,"AMPLITUDE"]
    amplitudes <- gain*amplitudes
    ranges <- reference$ranges[variables,]
    palette <- reference$palette

    # Restore attribute for binary variables.
    binary <- attr(results$planes, "binary")
    attr(comps, "binary") <- intersect(binary, variables)

    # Set colors and labels.
    colrs <- nroColorize(values=comps, amplitudes=amplitudes,
                         ranges=ranges, palette=palette)
    labls <- nroLabel(topology=topology, values=comps)

    # Launch a detached window.
    if(detach != "FALSE") {
        if(detach == "TRUE") grDevices::dev.new()
        if(detach == "aqua") {
	    if(capabilities("aqua")) grDevices::quartz()
	    else warning("Quartz display server not available.")
	}
        if(detach == "X11") {
	    if(capabilities("X11")) grDevices::x11()
	    else warning("X11 display server not available.")
	}
    }

    # Interactive subgrouping.
    try(topology <- nroPlot(topology=topology, colors=colrs, labels=labls,
        interactive=TRUE, clear=(detach == "FALSE")), silent=TRUE)

    # Convert to data frame.
    topology <- as.data.frame(topology, stringsAsFactors=FALSE)

    # Print report.
    t <- table(topology$REGION)
    if(sum(names(t) == "not_selected") < 1) {
        nsubs <- length(t)
        cat("\n", nsubs, " subgroups selected\n", sep="")
    } else {
        nsubs <- (length(t) - 1)
        cat("\n", nsubs, " + 1 subgroups selected\n", sep="")
    }
    return(topology)
}
