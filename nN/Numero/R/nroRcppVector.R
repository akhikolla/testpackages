nroRcppVector <- function(
    data,
    default,
    numeric=TRUE,
    empty=TRUE) {

    # Check if empty is allowed.
    if(length(data) < 1) {
        if(!empty) stop("Input is empty.")
        if(numeric) return(as.numeric(default))
        return(as.character(default))
    }

    # Check if input is a vector.
    if(!is.atomic(data) && !is.list(data)) {
	if(!empty) stop("Input is not a vector.")
        warning("Input is not a vector.")
        if(numeric) return(numeric())
        return(character())
    }

    # Copy element labels.
    labels <- names(data)

    # Set format.
    if(numeric) data <- as.numeric(data)
    else data <- as.character(data)

    # Return results.
    names(data) <- labels
    return(data)
}
