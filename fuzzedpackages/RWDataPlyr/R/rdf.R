
#' @export
print.rdf <- function(x, ...)
{
  mrm <- paste0("MRM: ", x$meta$name, "\n")
  desc <- paste0("Description: ", x$meta$description, "\n")
  n_trace <- paste0("Number of traces: ", x$meta$number_of_runs, "\n")
  obj <- paste0(paste0("   - ", names(x$runs[[1]]$objects)), collapse = "\n")
  
  cat(mrm, desc, n_trace, "Slots:\n", obj, sep = "")
  
  invisible(x)
}

#' @export
summary.rdf <- function(object, ...)
{
  cat("rdf contains", object$meta$number_of_runs, "traces of data,", 
      length(object$runs[[1]]$objects), "slots, and", 
      object$runs[[1]]$time_steps, "timesteps.")
}

#' Test if the object is an rdf
#' 
#' @param x An object
#' 
#' @return `TRUE` if the object inherits from the `rdf` class.
#' @export
is_rdf <- function(x) inherits(x, "rdf")

#' @rdname is_rdf
#' @export
is.rdf <- is_rdf
