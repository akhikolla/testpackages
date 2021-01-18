nroPostprocess <- function(
    data,
    mapping,
    reverse=FALSE,
    trim=FALSE) {

    # Nothing to do.
    if(length(mapping) < 1) return(data)

    # Check input.
    if(is.data.frame(mapping)) mapping <- attr(mapping, "mapping")
    if(is.matrix(mapping)) mapping <- attr(mapping, "mapping")

    # Set operation mode.
    if(reverse) {
        model.in <- mapping$output
        model.out <- mapping$input
    }
    else {
        model.in <- mapping$input
        model.out <- mapping$output
    }

    # Check model data.
    if(nrow(model.in) != nrow(model.out))
        stop("Incompatible model, size mismatch.")
    if(ncol(model.in) != ncol(model.out))
        stop("Incompatible model, size mismatch.")
    if(sum(rownames(model.in) != rownames(model.out)) > 0)
        stop("Incompatible model, row name conflict.")
    if(sum(colnames(model.in) != colnames(model.out)) > 0)
        stop("Incompatible model, column name conflict.")

    # Find variables.
    vars <- intersect(colnames(model.in), colnames(data))
    if(length(vars) < 1) {
        warning("No matching column names.")
        return(NULL)
    }

    # Prepare output.
    output <- data
    if(trim[[1]]) output <- NA*output

    # Process columns.
    for(vn in vars) {
        x <- model.in[,vn]
        y <- model.out[,vn]
        xout <- as.double(data[,vn])
        mask <- which(is.finite(x*y) & !duplicated(x))
        if(length(mask) < 3) next
        output[,vn] <- stats::approx(x=x[mask], y=y[mask],
	                   rule=2, xout=xout)$y
    }

    # Remove empty rows.
    if(trim[[1]]) {
        mu <- rowMeans(output, na.rm=TRUE)
        output <- output[which(is.finite(mu)),,drop=FALSE]
    }
    if(nrow(output) < 1) {
        warning("No usable rows.")
        return(NULL)
    }
    if(nrow(output) < nrow(data))
        warning("Unusable rows excluded.")

    # Remove empty columns.
    if(trim[[1]]) {
        mu <- colMeans(output, na.rm=TRUE)
        output <- output[,which(is.finite(mu)),drop=FALSE]
    }
    if(ncol(output) < 1) {
        warning("No usable columns.")
        return(NULL)
    }
    if(ncol(output) < ncol(data))
        warning("Unusable columns excluded.")

    # Return results.
    return(output)
}
