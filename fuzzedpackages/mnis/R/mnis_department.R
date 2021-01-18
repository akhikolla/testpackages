#' Request data on the holders of cabinet/shadow cabinet positions. Request specific departments by department ID (see \code{\link{mnis_reference}} for the \code{ref_department} function to retrieve departmental IDs).
#'
#' @param department_id The department look up. 0 returns the cabinet/shadow cabinet, -1 returns a list of all ministers. Defaults to 0.
#' @param bench Flag to return either Government or Opposition information. Defaults to 'Government'. The API is case sensitive on this parameter, so 'Government' or 'Opposition' will work, but 'government' and 'opposition' will not.
#' @param former Flag to include both current and former ministers/shadow ministers. Defaults to TRUE. If FALSE, only includes current ministers/shadow ministers.
#' @param tidy If TRUE, fixes the variable names in the tibble to remove non-alphanumeric characters and superfluous text, and convert to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @return A tibble with information on departments and ministers/shadow ministers.
#' @keywords mnis
#' @export
#' @seealso \code{\link{mnis_reference}}
#' @examples \dontrun{
#'
#' x <- mnis_department(department_id = 0, bench = 'Government', former=TRUE)
#'
#' }
#'


mnis_department <- function(department_id = 0, bench = "Government", former = TRUE, tidy = TRUE, tidy_style="snake_case") {

    if (former == TRUE) {

        query_former <- "former"

    } else {

      query_former <- "current"
    }

    department_id <- as.character(department_id)

    bench <- utils::URLencode(bench)

    baseurl <- "http://data.parliament.uk/membersdataplatform/services/mnis/Department/"

    query <- paste0(baseurl, department_id, "/", bench, "/", query_former, "/")

    got <- httr::GET(query, httr::accept_json())

    if (httr::http_type(got) != "application/json") {
        stop("API did not return json", call. = FALSE)
    }

    got <- mnis::tidy_bom(got)

    got <- jsonlite::fromJSON(got, flatten = TRUE)

    x <- tibble::as_tibble(as.data.frame(got$Department$Posts))

    if (tidy == TRUE) {

        x <- mnis::mnis_tidy(x, tidy_style)

        x

    } else {

        x

    }

}
