
#' Returns all members who are able to sit in either house, or who are currently ineligible to sit. Members ineligible to sit include but are not necessarily limited to former MPs, members of the judiciary, who are recused from House of Lords duties.
#' @param eligible If the member is currently eligible to sit. Accepts TRUE or FALSE. Defaults to TRUE.
#' @param house The house to which the member belongs. Accepts one of 'all', 'lords' and 'commons', defaults to 'all'. This parameter is not case sensitive, so 'commons', 'Commons' and 'cOmMOnS' will all return the same data.
#' @param party The party to which a member belongs. Defaults to NULL. The party must be fully spelled out (e.g. 'green party'), the API does not accept searches on this parameter. For a tibble of parties, see \code{\link{ref_parties}}. This parameter is not case sensititive.
#' @param tidy Fix the variable names in the tibble to remove special characters and superfluous text, and converts the variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case"
#' @keywords mnis
#' @export
#' @examples \dontrun{
#' x <- mnis_eligible(eligible=FALSE, house='all', party='labour')
#'
#' x <- mnis_eligible(eligible=TRUE, house='all', party='green party')
#'
#' x <- mnis_eligible(house="commons")
#'
#' }


mnis_eligible <- function(eligible = TRUE, house = "all", party = NULL, tidy = TRUE, tidy_style="snake_case") {

    if (is.na(pmatch(house, c("all", "lords", "commons"))))
        stop("Please select one of 'all', 'lords' or 'commons' for the parameter 'house'")

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/members/query/iseligible="

    house <- as.character(house)
    house <- tolower(house)

    if (house == "lords") {
        house <- "|house=lords"
    } else if (house == "commons") {
        house <- "|house=commons"
    } else if (house == "all") {
        house <- "|house=all"
    }

    if (is.null(party) == FALSE) {
        party <- utils::URLencode(party)
        party <- paste0("|party=", party)
    }

    query <- paste0(baseurl, eligible, house, party)

    got <- httr::GET(query, httr::accept_json(),encoding = "UTF-8")

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- tibble::as_tibble(got$Members$Member)

    if (tidy == TRUE) {

        x <- mnis::mnis_tidy(x, tidy_style)

        if(.Platform$OS.type=="windows"){

          x$member_from <- stringi::stri_trans_general(x$member_from, "latin-ascii")

          x$member_from <- gsub("Ynys MA\U00B4n", "Ynys M\U00F4n", x$member_from)

        }

        x

    } else {

        x

    }

}
