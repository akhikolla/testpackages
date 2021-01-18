
#' Function returns all members who took their seats in the house between two given dates.
#' @param start_date The start date of the search, Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Defaults to '1900-01-01'.
#' @param end_date The end date of the search. Accepts character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with \code{as.Date()}. Defaults to the current date.
#' @param house The house to which the member belongs. Accepts one of 'all', 'lords' and 'commons', defaults to 'all'.
#' @param party The party to which a member belongs. Defaults to NULL.
#' @param eligible If the member is currently eligible to sit. Accepts one of 'all', 'current', 'former', defaults to 'all'.
#' @param tidy Fix the variable names in the tibble to remove non-alphanumeric characters and superfluous text, and convert variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @keywords mnis
#' @export
#' @examples \dontrun{
#'
#' x <- mnis_joined_between(start_date = '2015-01-01', end_date ='2017-01-01', party='labour')
#'
#' }


mnis_joined_between <- function(start_date = "1900-01-01", end_date = Sys.Date(), house = "all", party = NULL, eligible = "all", tidy = TRUE, tidy_style="snake_case") {

    ## Making sure house works
    house <- as.character(house)

    house <- tolower(house)

    if (is.na(pmatch(house, c("all", "lords", "commons"))))
        stop("Please select one of 'all', 'lords' or 'commons' for the parameter 'house'")

    if (is.na(pmatch(eligible, c("all", "current", "former"))))
        stop("Please select one of 'all', 'current' or 'former' for the parameter 'eligible'")


    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/members/query/joinedbetween="

    start_date <- as.Date(start_date)

    end_date <- as.Date(end_date)

    if (is.null(party) == FALSE) {
        party <- utils::URLencode(party)
    }

    eligible <- as.character(eligible)

    if (house == "lords") {
        house <- "|house=lords"
    } else if (house == "commons") {
        house <- "|house=commons"
    } else if (house == "all") {
        house <- "|house=all"
    }

    if (is.null(party) == FALSE) {
        party <- paste0("|party*", party)
    }

    if (eligible == "all") {
        eligible <- NULL
    } else if (eligible == "current") {
        eligible <- "|iseligible=TRUE"
    } else if (eligible == "former") {
        eligible <- "|iseligible=FALSE"
    }

    query <- paste0(baseurl, start_date, "and", end_date, house, party, eligible)

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }
    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- tibble::as_tibble(got$Members$Member)

    if (tidy == TRUE) {

        x <- mnis::mnis_tidy(x, tidy_style)

        x

    } else {

        x

    }

}
