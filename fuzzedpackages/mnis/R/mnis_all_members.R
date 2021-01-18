
#' Returns a tibble with information on al members of both houses or a given house.
#'
#' @param house The house to which the member belongs. Accepts one of 'all', 'lords' and 'commons', defaults to 'all'. This parameter is not case sensitive, so using 'commons', 'Commons' or 'cOmMOnS' will result in the same data being returned.
#' @param party The party to which a member belongs. Defaults to NULL, in which case all members are returned, subject to other parameters. The party names are not case sensitive, but must be the complete string of the party name, searching and wildcard options are not accepted by the API, e.g. 'green party'.
#' @param tidy Fix the variable names in the tibble to remove special characters and superfluous text, and converts the variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case"
#' @keywords mnis
#' @return A tibble with information on all members of the House of Commons and/or the House of Lords that meet the criteria included in the function parameters.
#' @export
#'
#' @examples \dontrun{
#' x <- mnis_all_members(house = 'all', party = NULL, tidy = TRUE, tidy_style="snake_case")
#' }
#'

mnis_all_members <- function(house = "all", party = NULL, tidy = TRUE, tidy_style="snake_case") {

  house <- tolower(house)

  if (is.na(pmatch(house, c("all", "lords", "commons"))))
    stop("Please select one of 'all', 'lords' or 'commons' for the parameter 'house'")

  baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/members/query/Membership=all"

  if (is.null(party) == FALSE)
    party <- utils::URLencode(party)

  if (house == "lords") {
    house <- "|house=lords"
  } else if (house == "commons") {
    house <- "|house=commons"
  } else if (house == "all") {
    house <- ""
  }

  if (is.null(party) == FALSE) {
    party <- paste0("|party*", party)
  }

  message("Connecting to API")

  query <- paste0(baseurl, house, party, "/HouseMemberships/")

  got <- httr::GET(query, httr::accept_json())

  if (httr::http_type(got) != "application/json") {
    stop("API did not return json", call. = FALSE)
  }

  got <- mnis::tidy_bom(got)

  got <- jsonlite::fromJSON(got, flatten = TRUE)

  x <- got$Members$Member

  x <- tibble::as_tibble(x)

  if (tidy == TRUE) {

    x <- mnis::mnis_tidy(x, tidy_style)

    x

  } else {

    x

  }

}
