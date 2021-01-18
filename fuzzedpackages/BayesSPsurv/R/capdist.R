#' Gleditsch and Ward Distance data
#'
#' Dyadic dataset extracted from \href{https://journals.sagepub.com/doi/10.1177/0022343301038006006}{Gleditsch and Ward (2001)}.
#' The dataset contains information on the distace between capital cities among independent nation-states.
#'
#'
#'\describe{
#'    \item{numa}{COW code -- country A.}
#'    \item{ida}{Three letter ISO code -- country A.}
#'    \item{numb}{COW code -- country B.}
#'    \item{idb}{Three letter ISO code -- country B.}
#'    \item{kmdist}{Distance between capital cities in the kilometers.}
#'    \item{midist}{Minimal distance between capital cities in the kilometers.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name capdist
#' @usage data(capdist)
#' @format A data frame with 41006 rows and 6 variables
#' @source Gleditsch, Kristian S., and Michael D. Ward. (2001). "Measuring Space: A Minimum-Distance Database and Applications to International Studies." Journal of Conflict Resolution 38(6): 739-758.
NULL
