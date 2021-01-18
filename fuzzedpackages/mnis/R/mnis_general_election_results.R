
#' Returns an object containing list with details of the search parameter and a tibble with election results. Accepts queries on location type and name, and the start and end date to return general elections between. The API does not contain data for Norther Ireland.
#' @param location_type The type of area to return information for. Accepts 'Country', 'Region', 'County', and 'Constituency'. Defaults to 'Country'.
#' @param location_name The location to return data for. It can be the name of any Country, Region, County or Constituency. Defaults to 'Great Britain'.
#' @param start_date Start date of search. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Defaults to '1900-01-01' if no date is selected.
#' @param end_date End date of search. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Defaults to current date if no date is selected.
#' @param tidy Fix the variable names in the tibble to remove special characters and superfluous text, and converts the variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case"
#' @return Returns a list with details of the search parameter and a tibble with election results.
#' @keywords mnis
#' @export
#' @seealso \code{\link{mnis_reference}}
#' @examples \dontrun{
#' x <- mnis_general_election_results(location_type = 'Country', location_name = 'England',
#'                                  start_date = '2010-01-01', end_date = '2016-01-01')
#' }

mnis_general_election_results <- function(location_type = "Country", location_name = "Great Britain", start_date = "1900-01-01",
    end_date = Sys.Date(), tidy = TRUE, tidy_style="snake_case") {

    location_type <- utils::URLencode(location_type)

    location_name <- utils::URLencode(location_name)

    start_date <- as.Date(start_date)

    end_date <- as.Date(end_date)

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/GeneralElectionResults/"

    query <- paste0(baseurl, location_type, "/", location_name, "/", start_date, "/", end_date, "/")

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- got$ElectionResults

    x$ElectionResult <- as.tibble(x$ElectionResult)

    if (tidy == TRUE) {

        names(x)[names(x)=="LocationInfo"] <- "location_info"

        names(x)[names(x)=="ElectionResult"] <- "election_result"

        x$election_result <- mnis::mnis_tidy(x$election_result, tidy_style)

        x

    } else {

        x

    }

}
