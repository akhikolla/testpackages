#' Create a formula for renewalCount
#'
#' @param response the formula for the "main" parameter. It also specifies the
#'     response variable.
#' @param ... additional arguments for the ancilliary parameters.
#'
#' @return a Formula object suitable for argument \code{formula} of
#'     \code{renewalCount()}.
#' @export
CountrFormula <- function(response, ...){ # TODO: argument "dist"?
    dots <- list(...)
    wrk <- list(response, ...)

    res <- as.Formula(response, ...)
    if(length(dots) == 0)
        return(res)
    nams <- names(wrk)
    if(!is.null(nams)){
        dd <- rep(".", length(nams))
        ind <- which(nams != "")
        if(length(ind) > 0){
            if(length(ind) < length(dd) - 1){
                stop("All arguments in '...' must be named or all unnamed")
            }
            dd[ind] <- nams[ind]
        }
        dd <- as.formula(paste0(paste0(dd, collapse = " | "),  " ~ ." ))
        res <- update(res, dd)
    }
    res
}

