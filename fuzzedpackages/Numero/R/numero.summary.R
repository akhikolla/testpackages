numero.summary <- function(
    results,
    topology,
    data=NULL,
    capacity=10) {

    # Start processing.
    stamp <- date()
    cat("\n*** numero.summary ***\n", stamp, "\n", sep="")
    
    # Default dataset.
    if(is.null(data)) data <- results$data

    # Match results with original data.
    cat("\nResources:\n")
    layout <- results$layout
    pos <- match(rownames(layout), rownames(data))
    rows <- which(pos > 0)
    data <- data[pos[rows],]
    layout <- layout[rows,]
    cat(length(rows), " data points matched with layout\n", sep="")
    if(nrow(data) < 10) {
        cat("less than ten usable data points\n")
        return(NULL)
    }
    cat(ncol(data), " data columns\n", sep="")

    # Convert topology to data frame.
    topology <- data.frame(topology, stringsAsFactors=FALSE)
    if(is.null(topology$REGION)) {
        cat("no regions defined\n")
        return(NULL)
    }

    # Check subgroup capacity.
    nsubs <- length(table(topology$REGION))
    if(nsubs < 2) {
        cat("less than two subgroups\n")
        return(NULL)
    }
    if(length(t) > capacity) {
        cat("subgroup capacity exceeded\n")
        return(NULL)
    }

    # Check labels.
    if(is.null(topology$REGION.label)) {
        labls <- as.factor(topology$REGION)
        topology$REGION.label <- as.integer(labls)
	warning("Region labels set to defaults.")
    }

    # Estimate subgroup statistics.
    cat("\nComparisons:\n")
    suppressWarnings(
        output <- nroSummary(data=data, districts=layout$BMC,
                             regions=topology, capacity=capacity))
    if(length(output) < 1) {
        cat("no usable columns\n")
        return(output)
    }

    # Add region information to layout.
    layout$REGION <- attr(output, "regions")
    layout$REGION.label <- attr(output, "labels")
    attr(output, "layout") <- layout[,c("BMC","REGION","REGION.label")]
    attr(output, "regions") <- NULL
    attr(output, "labels") <- NULL

    # Find variables that had usable data.
    pvals <- output[, c("P.chisq", "P.t", "P.anova")]
    success <- which(rowMeans(pvals, na.rm=TRUE) >= 0)
    binary <- which(output$TYPE == "binary")
    categ <- which(output$TYPE == "categ")
    real <- which(output$TYPE == "real")
    binary <- unique(output$VARIABLE[intersect(success, binary)])
    categ <- unique(output$VARIABLE[intersect(success, categ)])
    real <- unique(output$VARIABLE[intersect(success, real)])
    cat(length(binary), " binary columns\n", sep="")
    cat(length(categ), " categorical columns\n", sep="")
    cat(length(real), " continuous columns\n", sep="")

    # Unusable variables.
    nskip <- (ncol(data) - length(binary) - length(categ) - length(real))
    cat(nskip, " unusable columns\n", sep="")
    return(output)
}
