nroDestratify <- function(
    data,
    labels) {

    # Convert input to numeric matrix.
    data <- nroRcppMatrix(data, trim=FALSE)

    # Check if anything to do.
    binary <- attr(data, "binary")
    numerics <- setdiff(colnames(data), binary)
    if(length(numerics) < ncol(data))
        warning("Binary or non-numeric columns.")
    if(length(numerics) < 1) {
        warning("No usable columns.")
        return(NULL)
    }

    # Check that inputs are compatible.
    labels <- nroRcppVector(labels, default=NULL, numeric=FALSE)
    grp <- as.integer(as.factor(labels))  
    if(nrow(data) != length(grp)) stop("Incompatible inputs.")

    # Check batch size.
    if(length(levels(grp)) > 0.2*length(grp))
        stop("Average batch size is less than five.")

    # Process only numeric columns.
    output <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
    rownames(output) <- rownames(data)
    colnames(output) <- colnames(data)
    data <- data[,numerics,drop=FALSE]

    # Remove differences in batch-specific distributions.
    res <- .Call("nro_destratify",
                 as.matrix(data),
                 as.integer(grp),
                 PACKAGE="Numero")  
    if(is.character(res)) stop(res)

    # Update values.
    names(res) <- colnames(data)
    for(vn in numerics)
        output[,vn] <- res[[vn]]

    # Determine columns that were not fully successful. 
    incomplete <- setdiff(colnames(output), numerics)
    for(vn in numerics) {
       xbits <- is.finite(data[,vn])
       ybits <- is.finite(output[,vn])
       if(sum(xbits) == sum(ybits)) next
       incomplete <- c(incomplete, vn)
    }
    if(length(incomplete) > 0)
        warning("Some values could not be processed.")
   
    # Update dataset.
    attr(output, "incomplete") <- incomplete
    return(output)
}
