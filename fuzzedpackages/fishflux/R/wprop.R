#' A function to find the ratio of dry weight and wet weight of fish in local database
#'
#' This function searches the ratio of dry weight and wet weight of fish
#' on the family level. If the family is not available, an average is returned.
#'
#' Returns a dataframe with the weight ratio (mdw) and it's sd (mdw_sd).
#'
#' @param family family
#' 
#' @keywords fish weight proportion
#' 
#' @importFrom stats median quantile sd
#'
#' @examples
#' library(fishflux)
#' wprop(family="Scaridae")
#' 
#' @export
wprop <- function(family) {
  wprop <- fishflux::weight_prop
  ww <- wprop[wprop$Family == family, "weight_prop"]
  ww_sd <- wprop[wprop$Family == family, "weight_prop_sd"]
  if (length(ww) == 0) {
    ww <- mean(wprop$weight_prop)
    ww_sd <- sd(wprop$weight_prop)
    warning("family not in database, average used")
  }
 data.frame(mdw = ww, mdw_sd = ww_sd)
}
