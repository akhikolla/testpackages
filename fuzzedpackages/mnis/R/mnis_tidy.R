
#' Functions to tidy up the variable names returned from the API, and turn dates and datetimes to POSIXct.
#' @param df The tibble to tidy.
#' @param tidy_style The style to tidy the tibble with.
#' @export
#' @rdname mnis_tidy
mnis_tidy <- function(df, tidy_style="snake_case") {

    df <- mnis::date_tidy(df)

    df <- member_tidy(df)

    names(df) <- gsub("@", "", names(df))

    names(df) <- gsub("#", "", names(df))

    names(df) <- gsub("\\.\\.", "\\.", names(df))

    names(df) <- gsub("^BasicDetails\\.", "", names(df))

    names(df) <- gsub("^BiographyEntries\\.", "", names(df))

    names(df) <- gsub("^Committees\\.", "", names(df))

    names(df) <- gsub("^Addresses\\.", "", names(df))

    names(df) <- gsub("^Constituencies\\.Constituency\\.", "Constituency.", names(df))

    names(df) <- gsub("^ElectionsContested\\.", "", names(df))

    names(df) <- gsub("^Experiences\\.", "", names(df))

    names(df) <- gsub("^GovernmentPosts\\.", "", names(df))

    names(df) <- gsub("^Honours\\.", "", names(df))

    names(df) <- gsub("^HouseMemberships\\.", "", names(df))

    names(df) <- gsub("^Statuses\\.", "", names(df))

    names(df) <- gsub("^Staff\\.", "", names(df))

    names(df) <- gsub("^Interests\\.Category\\.Interest\\.", "Interest\\.", names(df))

    names(df) <- gsub("^Interests\\.Category\\.", "Interest\\.", names(df))

    names(df) <- gsub("^MaidenSpeeches\\.", "", names(df))

    names(df) <- gsub("^OppositionPosts\\.", "", names(df))

    names(df) <- gsub("^Parties\\.", "", names(df))

    names(df) <- gsub("^PreferredNames\\.", "", names(df))

    names(df) <- gsub("^ParliamentaryPosts\\.", "", names(df))

    names(df) <- gsub("^OtherParliaments\\.", "", names(df))

    names(df) <- gsub("^ParliamentaryPosts\\.", "", names(df))

    names(df) <- gsub("^Post.PostHolders.PostHolder.Member", "PostHolder", names(df))

    names(df) <- gsub("^Post\\.PostHolders\\.", "", names(df))

    names(df) <- gsub("xsi:nil", "nil", names(df))

    names(df) <- gsub("xmlns:xsi", "label", names(df))

    names(df) <- gsub("\\.", "_", names(df))

    names(df) <- gsub("([[:lower:]])([[:upper:]])", "\\1_\\2", names(df))

    names(df) <- gsub("__", "_", names(df))

    names(df) <- tolower(names(df))

    names(df)[names(df) == "df_house"] <- "house"

    names(df) <- gsub("_xsi:nil", "_nil", names(df))

    names(df) <- gsub("_xmlns:xsi", "", names(df))

    if(tidy_style=="camelCase") {

      names(df) <- gsub("(^|[^[:alnum:]])([[:alnum:]])", "\\U\\2", names(df), perl = TRUE)

      substr( names(df), 1, 1) <- tolower(substr(names(df), 1, 1))

    } else if (tidy_style=="period.case") {

      names(df) <- gsub("_", ".", names(df))

    }

    df

}

#' member_tidy
#' @export
#' @rdname mnis_tidy
member_tidy <- function(df){

  names(df) <- gsub("^Members\\.Member\\.", "", names(df))

  df

}


#' ref_tidy
#' @export
#' @rdname mnis_tidy
ref_tidy <- function(df, tidy_style) {

    names(df) <- gsub(".*\\.", "", names(df))

    names(df) <- gsub("([[:lower:]])([[:upper:]])", "\\1_\\2", names(df))

    names(df) <- tolower(names(df))

    if(tidy_style=="camelCase") {

      names(df) <- gsub("(^|[^[:alnum:]])([[:alnum:]])", "\\U\\2", names(df), perl = TRUE)

      substr( names(df), 1, 1) <- tolower(substr(names(df), 1, 1))

    } else if (tidy_style=="period.case") {

      names(df) <- gsub("_", ".", names(df))

    }

    df

}

#' constituency_results_tidy
#' @param results The tibble to tidy
#' @param details The list to tidy
#' @export
constituency_results_tidy <- function(results, details) {

    names(results) <- gsub("Candidates\\.Candidate\\.", "", names(results))

    names(results) <- gsub("\\.", "", names(results))

    names(results) <- gsub("([[:lower:]])([[:upper:]])", "\\1_\\2", names(results))

    names(details) <- gsub("([[:lower:]])([[:upper:]])", "\\1_\\2", names(details))

    names(results) <- tolower(names(results))

    names(details) <- tolower(names(details))

    y <- c(list(results = results), list(details = details))

    y

}


#' A function to strip Byte Order Marks (BOM) out of JSON data returned from the API.
#'
#' @param df The GET returned from call to API.
#' @export
tidy_bom <- function(df) {

    got <- as.character(df)

    got <- mnis_bom(got)

    got

}


#' A function that makes date variables returned from the API datable, ie by converting them to POSIXct. Does the same thing for datetimes.
#'
#' @param df The tibble with the undateable dates.
#' @export
date_tidy <- function(df) {

  date_vars <- grepl('date', colnames(df), ignore.case=TRUE)

  df[date_vars] <- lapply(df[date_vars], function(y) gsub("T", " ", y))

  df[date_vars] <- lapply(df[date_vars], as.POSIXct, format = "%Y-%m-%d %H:%M:%S")

  df

}
