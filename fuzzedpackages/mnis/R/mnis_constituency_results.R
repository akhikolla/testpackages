
#' Returns a single object containing a list with details of the constituency and a tibble with election results. For constituency IDs, see \code{\link{ref_constituencies}}.
#' @param constituency_id The ID of the constituency to return the data for. This parameter cannot be empty.
#' @param election_id The ID of the election to return the data for. Defaults to 0, which calls the most recent result, either the result of the last general election, or the result of the last byelection held since that election.
#' @param tidy If TRUE, fixes the variable names in the tibble to remove non-alphanumeric characters and superfluous text, and convert to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @return A list with details of the constituency, labelled 'details' and a tibble with election results, labelled 'results'. The list and tibble are stored in a single object.
#' @keywords mnis
#' @export
#' @seealso \code{\link{mnis_reference}}
#' @examples \dontrun{
#'
#' x <- mnis_constituency_results(constituency_id = 3709, election_id = 0)
#'
#' }

mnis_constituency_results <- function(constituency_id = NULL, election_id = 0, tidy = TRUE, tidy_style="snake_case") {

    if (missing(constituency_id)) {
        stop("'constituency_id' cannot be empty", call. = FALSE)
    }

    constituency_id <- as.character(constituency_id)

    election_id <- as.character(election_id)

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/ConstituencyResults/"

    query <- paste0(baseurl, constituency_id, "/", election_id, "/")

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    details <- got$Constituency$Details

    results <- tibble::as_tibble(got$Constituency$Results)

    y <- list()

    if (tidy == TRUE) {

        y <- mnis::constituency_results_tidy(results, details)

        y

    } else {

        y <- c(list(results = results), list(details = details))

        y

    }

}
