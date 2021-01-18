
#' A series of basic function for the API lookup. Each function accepts a member's ID and returns information; if no ID is given basic information on all members of both houses is returned.
#' @param ID The member ID value. If empty, function calls \code{\link{mnis_all_members}} and returns basic information on all members of both houses.
#' @param ref_dods Request based on the DODS membership ID scheme. Defaults to FALSE, where it requests data based on the default membership ID scheme.
#' @param tidy If TRUE, fixes the variable names in the tibble to remove non-alphanumeric characters and superfluous text, and convert to a consistent style. Defaults to TRUE.
#' @param tidy_style The style to convert variable names to, if tidy=TRUE. Accepts one of "snake_case", "camelCase" and "period.case". Defaults to "snake_case".
#' @keywords mnis
#' @return A tibble with the data corresponding to the particular function called.
#' @examples \dontrun{
#'
#' x <- mnis_basic_details(172)
#'
#' }
#' @export
#' @rdname mnis_additional
#' @seealso \code{\link{mnis_full_biog}} \code{\link{mnis_extra}}
#'
#' @export
#' @rdname mnis_additional
#' @examples \dontrun{
#'
#' x <- mnis_additional()
#'
#' }

mnis_additional <- function() {

    x <- c("mnis_full_biog()", "mnis_basic_details()", "mnis_biography_entries()", "mnis_committees()", "mnis_addresses()", "mnis_constituencies()", "mnis_elections_contested()", "mnis_experiences()", "mnis_government_posts()", "mnis_honours()", "mnis_house_memberships()", "mnis_statuses()", "mnis_staff()", "mnis_interests()", "mnis_known_as()", "mnis_maiden_speeches()", "mnis_opposition_posts()", "mnis_other_parliaments()", "mnis_parliamentary_posts()", "mnis_parties()", "mnis_preferred_names()")
    message("All Available Additional Information Functions:")

    print(x)

}

#' @export
#' @rdname mnis_additional
mnis_basic_details <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/BasicDetails")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional

mnis_biography_entries <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/BiographyEntries")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_committees <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Committees")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_addresses <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Addresses")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_constituencies <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Constituencies")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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
#' @export
#' @rdname mnis_additional
mnis_elections_contested <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/ElectionsContested")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_experiences <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Experiences")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_government_posts <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/GovernmentPosts")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_honours <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Honours")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_house_memberships <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/HouseMemberships")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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
#' @export
#' @rdname mnis_additional

mnis_statuses <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Statuses")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional

mnis_staff <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Staff")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_interests <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Interests")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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
#' @export
#' @rdname mnis_additional
mnis_known_as <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/KnownAs")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional
mnis_maiden_speeches <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/MaidenSpeeches")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional

mnis_opposition_posts <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/OppositionPosts")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional

mnis_other_parliaments <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/OtherParliaments")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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
#' @export
#' @rdname mnis_additional
mnis_parliamentary_posts <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/ParliamentaryPosts")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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

#' @export
#' @rdname mnis_additional

mnis_parties <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/Parties")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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
#' @export
#' @rdname mnis_additional

mnis_preferred_names <- function(ID = NULL, ref_dods = FALSE, tidy = TRUE, tidy_style="snake_case") {

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

        query <- paste0(baseurl, ID_Type, ID, "/PreferredNames")

        got <- httr::GET(query, httr::accept_json())

        if (httr::http_type(got) != "application/json") {
            stop("API did not return json", call. = FALSE)
        }

        got <- mnis::tidy_bom(got)

        got <- jsonlite::fromJSON(got, flatten = TRUE)

        # got <- jsonlite::fromJSON(httr::content(got, 'text', encoding = 'bytes'), flatten = TRUE)

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
