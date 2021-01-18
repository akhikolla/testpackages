#' List of countries available in the Covid-19 database
#'
#' @description This function provides a list of 183 country names inside the Covid-19 database.
#' The \code{country_name} argument inside the \code{\link{load_covid}} function should be spelled exactly
#' as they are listed here. The user must be online for this function to work.
#'
#' @source \url{https://github.com/CSSEGISandData/COVID-19}
#'
#' @seealso \code{\link{load_covid}}
#'
#' @importFrom curl has_internet
#' @importFrom utils read.csv
#' @export

country_list <- function(){

  if(curl::has_internet() == F) stop("The user must be online to use this function")

  baseURL = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series"
  covid <- utils::read.csv(file.path(baseURL,"time_series_covid19_confirmed_global.csv"), check.names=FALSE, stringsAsFactors=FALSE)
  covid <- dplyr::rename(covid,country = 'Country/Region')
  covid <- dplyr::filter(covid,!(covid$country %in% c("Diamond Princess","Holy See","Taiwan*","MS Zaandam","Western Sahara")))

  country_list <- unique(covid$country)
  sort(country_list)

}
