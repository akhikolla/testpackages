
#' Requests all available biographical information for a given member, and returns it in the form of a tibble.
#' @param ID The ID number of the member, using the default MNIS scheme. If \code{ref_dods} is TRUE, accepts the Dods monitoring scheme instead. If left empty, returns the same data as \code{\link{mnis_all_members}} with default parameters.
#' @param ref_dods Request based on the Dods monitoring member ID scheme. Defaults to FALSE. If FALSE, requests using the default MNIS identification scheme.
#' @param tidy Fix the variable names in the tibble to remove non-alphanumeric characters and superfluous text, and convert variable names to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @keywords mnis
#' @export
#' @examples \dontrun{
#'
#' x <- mnis_full_biog(172)
#'
#' }
#' @seealso \code{\link{mnis_basic_details}} \code{\link{mnis_additional}} \code{\link{mnis_extra}}

mnis_full_biog <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

    if (missing(ID)) {

        x <- mnis_all_members()

    } else {

        ID <- as.character(ID)

        if (ref_dods == TRUE) {
            ID_Type <- "refDods="
        } else {
            ID_Type <- "id="
        }

        baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/members/query/"

        query <- paste0(baseurl, ID_Type, ID, "/FullBiog")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        dl <- data.frame(ID = rep(names(got), sapply(got, length)), Obs = unlist(got))

        x <- t(dl)

        x <- as.data.frame(x)

        x <- x[rownames(x) != "ID", ]

        x <- tibble::as_tibble(x)

    }

    if (tidy == TRUE) {

      x <- mnis::mnis_tidy(x, tidy_style)

      x

    } else {

        x

    }

}
