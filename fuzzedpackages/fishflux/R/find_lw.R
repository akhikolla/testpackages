#' A function to find length-weight relationship parameters a and b
#'
#' A function to find estimates length-weight relationship parameters available on fishbase. It returns a list of means and standard deviations of a and b obtained from:
#'*Froese, R., J. Thorson and R.B. Reyes Jr., 2013. A Bayesian approach for estimating length-weight relationships in fishes. J. Appl. Ichthyol. (2013):1-7.*
#'Please cite Froese et al. (2013), when using these values.
#'The default mirror for fishbase is set to "de", please change this if needed for your location
#'
#' @param sp A charachter value containing the species name
#' @param mirror Mirror for fishbase (eg. "de", "org", "us", etc.) Default is "us".
#' 
#' @keywords fish l-w relationship fishbase bayesian
#' 
#' @importFrom httr GET
#' @importFrom httr message_for_status
#' @importFrom curl has_internet
#' 
#' @examples
#' \donttest{library(fishflux)
#' library(plyr)
#' # find length-weight relationship parameters for one species
#' find_lw("Lutjanus griseus")
#'
#' # find length-weight relationship parameters for multiple species and return in dataframe
#' ldply(lapply(c("Chlorurus spilurus","Zebrasoma scopas"), find_lw))}
#' 
#' @export
find_lw <- function(sp, mirror = "us") {

  check_name_fishbase(sp)

  sp_ <- gsub(" ", "-", sp)
  url <- paste("https://www.fishbase.", mirror, "/summary/", sp_, ".html", sep = "")
  
  try_GET <- function(x, ...) {
    tryCatch(
      httr::GET(url = x, httr::timeout(10), ...),
      error = function(e) conditionMessage(e),
      warning = function(w) conditionMessage(w)
    )
  }
  is_response <- function(x) {
    class(x) == "response"
  }
  
  # First check internet connection
  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  # Then try for timeout problems
  resp <- try_GET(url)
  if (!is_response(resp)) {
    message(resp)
    return(invisible(NULL))
  }
  # Then stop if status > 400
  if (httr::http_error(resp)) { 
    httr::message_for_status(resp)
    return(invisible(NULL))
  }
  
  # ready to read page
  page <-  readLines(url)

  l <- grep("Bayesian length-weight:", page)
  line <- page[l]
  line <- gsub("\t", "", line)
  line <- gsub("Bayesian length-weight: ", "", line)
  line <- gsub(", in cm Total Length, based on LWR estimates for this species (Ref. <A href='../references/FBRefSummary.php?ID=93245'>93245</A>).", "", line)
  ob <- strsplit(line, " ")
  ob <- ob[[1]]

  lwa_m <- as.numeric(gsub("a=", "", ob[1]))
  lwa_up <- as.numeric(gsub("),", "", ob[4]))
  lwa_sd <- (lwa_up - lwa_m) / 1.96

  lwb_m <- as.numeric(gsub("b=", "", ob[5]))
  lwb_up <- as.numeric(gsub("),", "", ob[8]))
  lwb_sd <- (lwb_up - lwb_m) / 1.96

  data.frame(species = sp, lwa_m = lwa_m,
             lwa_sd = lwa_sd, lwb_m = lwb_m,
             lwb_sd = lwb_sd)
}
