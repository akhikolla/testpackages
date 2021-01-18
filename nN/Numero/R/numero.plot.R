numero.plot <- function( 
    results,
    variables=NULL,
    topology=NULL,
    folder=NULL,
    prefix="figure",
    reference=NULL,
    subplot=NULL,
    gain=1.0,
    detach=FALSE,
    capacity=500,
    font=NULL) {

    # Start processing.
    cat("\n*** numero.plot ***\n", date(), "\n", sep="")

    # Default inputs.
    if(is.null(variables)) variables <- colnames(results$planes)
    if(is.null(topology)) topology <- results$map$topology
    if(is.null(reference)) reference <- results
    prefix <- as.character(prefix[[1]])
    detach <- as.character(detach[[1]])

    # Check capacity.
    cat("\nResources:\n")
    if(capacity < 2) {
        cat("capacity less than two\n", sep="")
        return(0)
    }

    # Select variables.
    comps <- results$planes
    stats <- results$statistics
    variables <- intersect(variables, colnames(comps))
    cat(length(variables), " columns included\n", sep="")
    if(length(variables) < 1) {
        cat("no usable variables\n", sep="")
        return(0)
    }

    # Check if too many variables.
    stats <- stats[variables,]
    comps <- comps[,variables,drop=FALSE]
    if(nrow(stats) > capacity) {
        cat("capacity exceeded, showing ", capacity, " plots.\n", sep="")
        comps <- comps[,1:capacity,drop=FALSE]
	stats <- stats[1:capacity,]
    }

    # Check if reference is usable.
    variables <- colnames(comps)
    rvars <- rownames(reference$statistics)
    if(sum(is.na(match(variables, rvars))) > 0) {
	warning("Incompatible reference.")
        variables <- intersect(rvars, variables)
        if(length(variables) < 1) return(0)
        comps <- comps[,variables]
        stats <- stats[variables,]
    }

    # Check if gain is usable.
    gain <- as.double(gain[[1]])
    if(!is.finite(gain)) {
        gain <- 1
        cat("unusable gain, reverted to 1\n", sep="")
    }
    if(gain <= 0.0) {
        gain <- 1
        cat("non-positive gain, reverted to 1\n", sep="")
    }

    # Check if subplot is usable.
    if(length(subplot) < 2) {
        if(is.null(folder)) subplot <- c(3,3)
	else subplot <- c(10,4)
    }
    subplot <- as.integer(subplot[1:2])
    if((subplot[1] < 1) || (subplot[2] < 1)) {
        cat("unusable subplot, reverting to automatic\n", sep="")
        if(length(folder) < 1) subplot <- c(3,3)
	else subplot <- c(10,4)
    }

    # Check if folder is accessible.
    if(length(folder) > 0) {
        if(!dir.exists(folder)) dir.create(folder)
	if(!dir.exists(folder)) {
	    cat("destination '", folder, "' not available\n", sep="")
	    folder <- NULL
	}
	if(!is.null(folder))
	    cat("destination folder '", folder, "'\n", sep="")
    }
    else {
        cat("destination folder not defined\n", sep="")
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

    # Extract attributes.
    contrast <- attr(colrs, "contrast")
    visible <- attr(labls, "visible")

    # Split into several figures.
    nfigs <- 0
    nsubs <- (subplot[1])*(subplot[2])
    nstats <- nrow(stats)
    while(TRUE) {

        # Select colorings.
        mask <- (nfigs*nsubs + 1:nsubs)
        mask <- mask[which(mask <= nstats)]
	if(length(mask) < 1) break

        # Print progress message.
        nfigs <- (nfigs + 1)
        cat("\nFigure ", nfigs, ":\n", sep="")
        cat(length(mask), " subplots\n", sep="")

        # Set file names.
        fn.svg <- NULL
        fn.html <- NULL
        if(length(folder) > 0) {
	    fn.svg <- sprintf("%s%02d.svg", prefix, nfigs)
	    fn.svg <- file.path(folder, fn.svg)
	    fn.html <- sprintf("%s%02d.html", prefix, nfigs)
	    fn.html <- file.path(folder, fn.html)
	    cat("file name '", fn.svg, "'\n", sep="")
        }

        # Make sure column names are preserved.
        colrs.masked <- as.matrix(colrs[,mask])
        labls.masked <- as.matrix(labls[,mask])
        comps.masked <- as.matrix(comps[,mask])
        if(length(mask) == 1) {
            cname <- colnames(colrs)
	    if(length(cname) < 1) cname <- mask
	    else cname <- cname[mask]
	    colnames(colrs.masked) <- cname
	    colnames(labls.masked) <- cname
	    colnames(comps.masked) <- cname
        }

        # Restore attributes.
        attr(colrs.masked, "contrast") <- as.matrix(contrast[,mask])
        attr(labls.masked, "visible") <- as.matrix(visible[,mask])

        # Launch a detached window.
        if((length(fn.svg) < 1) && (detach != "FALSE")) {
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

        # Plot colorings.
	if(is.null(fn.svg)) {
             nroPlot(topology=topology, colors=colrs.masked,
	         labels=labls.masked, subplot=subplot)
	     next
	}

        # Save colorings.
        nbyt <- nroPlot.save(file=fn.svg, topology=topology, font=font,
	    colors=colrs.masked, labels=labls.masked, subplot=subplot)
        cat(nbyt, " bytes saved in '", fn.svg, "'\n", sep="")
        nbyt <- nroPlot.save(file=fn.html, topology=topology, font=font,
	    colors=colrs.masked, labels=labls.masked, subplot=subplot)
        cat(nbyt, " bytes saved in '", fn.html, "'\n", sep="")
    }

    # Final report.
    cat("\nSummary:\n")
    if(length(folder) < 1) cat(nfigs, " figures\n", sep="")
    else cat(nfigs, " figures -> '", folder, "'\n", sep="")
    return(nfigs)
}
