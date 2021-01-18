#' Simple aggregation functions for monthly matrix data
#' 
#' A family of functions that take a matrix containing monthly data (months by 
#' traces) that has a "timespan" attribute, annualizes the data by summing, or
#' finding the minimum or maximum monthly values. Returns a years by traces 
#' matrix. Matrices returned by [rdf_get_slot()] have the timespan attribute 
#' added to them. 
#' 
#' @param rwslot A matrix (months by traces) such as that returned by 
#'   [rdf_get_slot()]. Function will error if the `rwslot` does not contain 
#'   "regular" monthly data, i.e., the data must start in January and end in
#'   December, or start in October and end in September (water year), and the 
#'   `rwslot` must have the timespan attribute. 
#'   
#' @param multFactor A factor the annual sum will be multiplied by.  Can be used 
#'   to convert from flow to volume, or to scale all results in another manor.
#'   
#' @return Other functions: Annual matrix (years x traces)
#'   
#' @seealso [rdf_get_slot()]
#' 
#' @examples
#' zz <- rdf_get_slot(keyRdf, 'Powell.Outflow')
#' 
#' # returns in original units, e.g., acre-ft
#' annualTotVal <- rwslot_annual_sum(zz)
#' 
#' # returns in scaled units, e.g., kaf
#' annualTotVal <- rwslot_annual_sum(zz, 0.001) 
#' 
#' @rdname rwslot_aggs
#' @export
rwslot_annual_sum <- function(rwslot, multFactor = 1) {
  # take each column, make it a matrix of years by
  check_rwslot(rwslot, as.character(sys.call(sys.parent()))[1L])
  
  res <- do.call(
    cbind, 
    lapply(
      seq_len(ncol(rwslot)),
      function(xx) apply(matrix(rwslot[,xx],ncol = 12, byrow = TRUE), 1, sum)
    )
  )
  colnames(res) <- colnames(rwslot)
  res * multFactor
}

#' @describeIn rwslot_aggs Deprecated version of `rwslot_annual_sum()`.
#' @export
sumMonth2Annual <- function(rwslot, multFactor = 1) 
{
  .Deprecated("rwslot_annual_sum()")
  rwslot_annual_sum(rwslot, multFactor)
}

#' @describeIn rwslot_aggs finds the minimum annual value 
#'   for all years and traces.
#'   
#' @examples
#' pe <- rdf_get_slot(keyRdf,'Mead.Pool Elevation')
#' peMax <- rwslot_annual_min(pe)
#' 
#' @export
rwslot_annual_min <- function(rwslot)
{
  check_rwslot(rwslot, as.character(sys.call(sys.parent()))[1L])
  apply(rwslot, 2, trace_min_ann)
}

#' @export
#' @describeIn rwslot_aggs Deprecated version of `rwslot_annual_min()`.
getMinAnnValue <- function(rwslot)
{
  .Deprecated("rwslot_annual_min()")
  rwslot_annual_min(rwslot)
}

trace_min_ann <- function(traceVal)
{
  tmp <- matrix(traceVal, ncol = 12, byrow = TRUE)
  apply(tmp, 1, min)
}

#' @describeIn rwslot_aggs finds the maximum annual value for all years and 
#'   traces.
#' 
#' @examples
#' pe <- rdf_get_slot(keyRdf,'Mead.Pool Elevation')
#' peMax <- rwslot_annual_max(pe)
#' 
#' @export
rwslot_annual_max <- function(rwslot)
{
  check_rwslot(rwslot, as.character(sys.call(sys.parent()))[1L])
  apply(rwslot, 2, trace_max_ann)
}

#' @export
#' @describeIn rwslot_aggs Deprecated version of `rwslot_annual_max()`.
getMaxAnnValue <- function(rwslot)
{
  .Deprecated("rwslot_annual_max()")
  rwslot_annual_max(rwslot)
}

trace_max_ann <- function(traceVal)
{
  tmp <- matrix(traceVal, ncol = 12, byrow = TRUE)
  
  apply(tmp, 1, max)
}

#' @describeIn rwslot_aggs calculates the flow-weighted average annual 
#' concentration (fwaac). Given mass and flow at the monthly basis, the 
#' flow-weighted average annual concentration is computed. `mass` and `flow` 
#' should be monthly data. `rwslot_fwaac()` expects flow to be in acre-ft/month 
#' and mass to be in tons; however, there are no checks to ensure this is true.
#' Return value will be in mg/L.
#' 
#' @param mass A matrix (months by traces), such as that returned by 
#'   [rdf_get_slot()], of mass in tons.
#' @param flow A matrix (months by traces), such as that returned by 
#'   [rdf_get_slot()], of flow in acre-ft/month.

#' @return `rwslot_fwaac()`: Annual matrix (years x traces). Units are mg/L.
#'   
#' @examples
#' flow <- rdf_get_slot(keyRdf,'Powell.Outflow')
#' # make up mass, since it's not stored in the example data
#' rr <- matrix(
#'   rnorm((nrow(flow) * ncol(flow)), mean = 1000, sd = 200), 
#'   nrow = nrow(flow), 
#'   ncol = ncol(flow)
#' )
#' mass <- flow / 1000000 * rr^2 - rr + 1500 
#' fwaac <- rwslot_fwaac(mass, flow) 
#' 
#' @export
rwslot_fwaac <- function(mass, flow)
{
  if (!identical(dim(mass), dim(flow))) {
    stop("In `rwslot_fwaac()`, the dimensions of `flow` and `mass` must match.")
  }
  
  check_rwslot(mass, as.character(sys.call(sys.parent()))[1L])
  check_rwslot(flow, as.character(sys.call(sys.parent()))[1L])
  
  nyear <- nrow(mass)/12
  
  vapply(
    seq_len(ncol(mass)), 
    function(x) trace_fwaac(mass[,x], flow[,x]), 
    FUN.VALUE = numeric(nyear)
  )
}

trace_fwaac <- function(mass, flow)
{
  if (!is_full_monthly(length(mass)) || !is_full_monthly(length(flow))) {
    stop('Data passed to `trace_fwaac()` is not divisible by 12')
  }
  # move into a years x months matrix
  mass <- matrix(mass, ncol = 12, byrow = TRUE)
  flow <- matrix(flow, ncol = 12, byrow = TRUE)
  
  mass.annAvg <- apply(mass, 1, sum)/12
  flow.annAvg <- apply(flow, 1, sum)/12 # now essentially a volume
  
  conc <- mass.annAvg/flow.annAvg*735.466642
  
  conc
}

#' is_full_monthly assumes that data is monthly if it is divisible by 12
#' @noRd
is_full_monthly <- function(x)
{
  x%%12 == 0
}

#' is_regular_timespan checks to see if it is "regular" monthly data, i.e., 
#' that it starts in January and ends in December, or starts in October and
#' ends in September
#' @param rwslot Matrix from [rdf_get_slot()]
#' @noRd
is_regular_timespan <- function(rwslot)
{
  ts_months <- ym_get_month_str(zoo::as.yearmon(attr(rwslot,"timespan")))
  
  all(ts_months == c("January", "December")) || 
    all(ts_months == c("October", "September"))
}

#' Checks that the rwslot timespan is regular monthly data, and erros if it's
#' not
#' @noRd
check_rwslot <- function(rwslot, func)
{
  if (!(is.matrix(rwslot) && "timespan" %in% names(attributes(rwslot))))
    stop(
      "`", func, "()`, expects a matrix with a 'timespan' attribute.",
      call. = FALSE
    )
  
  if (!is_regular_timespan(rwslot))
    stop(
      "`", func, "()`, expects a regular monthly timespan.\n",
      "I.e., it should start in January and end in December or start in\n",
      "October and end in September.",
      call. = FALSE
    )
  
  invisible(rwslot)
}
