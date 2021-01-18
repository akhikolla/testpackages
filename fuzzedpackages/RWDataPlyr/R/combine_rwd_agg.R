
#' Combine RiverWare data aggregators
#' 
#' Take a sequence of `rwd_agg` arguments (or vector, matrix, or data.frames) 
#' and combine by rows. If the objects are not `rwd_agg` objects they will be 
#' combined through the default `rbind()` method, and then verified that they
#' meet all constraints to be a valid `rwd_agg` object. `cbind()` will fail for
#' `rwd_agg` objects.
#' 
#' @inheritParams base::cbind
#' @export
#' @examples 
#' 
#' ra1 <- rwd_agg(data.frame(
#'   file = "KeySlots.rdf",
#'   slot = "Powell.Pool Elevation",
#'   period = "wy",
#'   summary = "min",
#'   eval = "<",
#'   t_s = 3550,
#'   variable = "powellLt3550",
#'   stringsAsFactors = FALSE
#' ))
#' 
#' ra2 <- read_rwd_agg(
#'   system.file(
#'     "extdata/rwd_agg_files/passing_aggs.csv", 
#'     package = "RWDataPlyr"
#'   )
#' )
#' 
#' rbind(ra1, ra2)
#' 
#' \dontrun{
#' # will fail because you cannot have repeating variable names
#' rbind(ra1, ra1)
#' 
#' # will also fail
#' cbind(ra1, ra2)
#' }
#' 
rbind.rwd_agg <- function(..., deparse.level = 1)
{
  # couldn't get UseNextMethod to work, but this works
  validate_rwd_agg(rbind.data.frame(..., stringsAsFactors = FALSE))
}

#' @inheritParams base::cbind
#' @export
#' @rdname rbind.rwd_agg
cbind.rwd_agg <- function(..., deparse.level = 1)
{
  stop("`rwd_agg` objects cannot be combined by columns.")
}
