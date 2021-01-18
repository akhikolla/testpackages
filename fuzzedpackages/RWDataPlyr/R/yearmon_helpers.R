
#' Get the water year from a year-month (yearmon) value
#' 
#' `ym_get_wateryear()` returns the water year (assumed to be October - 
#' September) from a [zoo::yearmon] object.
#' 
#' If the argument is not already a yearmon object, it will attempt to convert 
#' it to a [zoo::yearmon]. This may result in unexpected results. For example, 
#' the string `"12-1-1906"` can be converted to a [zoo::yearmon], however, it 
#' will not convert to `"Dec 1906"` as you might desire. It will convert to 
#' `"Jan 0012"` since it is not a format expected by [zoo::as.yearmon()]. 
#' Therefore, a warning is posted when the function attempts to convert to 
#' [zoo::yearmon], and it is safer to ensure `ym` is already a [zoo::yearmon]. 
#' 
#' @param ym An object of class [zoo::yearmon], or something that can be 
#'   successfully converted to [zoo::yearmon].
#'   
#' @return The water year as a numeric.
#' 
#' @examples 
#' ym_get_wateryear(zoo::as.yearmon(c("Dec 1906", "Oct 1945", "Jul 1955")))
#' ym_get_wateryear("2000-11")
#' 
#' @export

ym_get_wateryear <- function(ym)
{
  if (!methods::is(ym, "yearmon")) {
    warning("ym, is not a yearmon object. attempting to convert to yearmon...")
    ym <- zoo::as.yearmon(ym)
    if (is.na(ym)) 
      stop("could not convert ym to yearmon")
  }
  
  mm <- as.numeric(format(ym, '%m'))
  yy <- ym_get_year(ym)
  # if OND then increase year by one for water year, else leave it the same
  yy[mm > 9] <- yy[mm > 9] + 1
  
  yy
}

#' @export
#' @rdname ym_get_wateryear
getWYFromYearmon <- function(ym)
{
  .Deprecated("ym_get_wateryear()")
  ym_get_wateryear(ym)
}

#' Get the year as a numeric from a yearmon object
#' 
#' Could use lubridate::year(), but for now we are not depending on lubridate
#' @noRd
ym_get_year <- function(ym)
{
  if (!(class(ym) == "yearmon"))
    stop("ym in ym_get_year(ym) is not a yearmon object.")
  
  as.numeric(format(ym, "%Y"))
}

#' Get the full month name from a yearmon object
#' @noRd
ym_get_month_str <- function(ym)
{
  if (!(class(ym) == "yearmon"))
    stop("ym in ym_get_month_str(ym) is not a yearmon object.")
  
  format(ym, "%B")
}
