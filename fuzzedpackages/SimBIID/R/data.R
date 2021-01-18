#' @title Time series counts of ebola cases
#'
#' @description A dataset containing time series counts for the number
#'      of new individuals exhibiting clinical signs, and the number of
#'      new removals each day for the 1995 Ebola epidemic in the Democratic
#'      Republic of Congo
#'
#' @format A data frame with 192 rows and 3 variables:
#' \describe{
#'   \item{time}{days from 1st January 1995}
#'   \item{clin_signs}{number of new clinical cases at each day}
#'   \item{removals}{number of new removals at each day}
#' }
#' @source Khan AS et al. (1999) <doi:10.1086/514306>
"ebola"

#' @title Time series counts of smallpox cases
#'
#' @description A dataset containing time series counts for the number
#'      of new removals for the 1967 Abakaliki smallpox outbreak.
#'
#' @format A data frame with 23 rows and 2 variables:
#' \describe{
#'   \item{time}{days from initial observed removal}
#'   \item{removals}{number of new removals in (time - 1, time)}
#' }
#' @source Thompson D and Foege W (1968) <https://apps.who.int/iris/bitstream/handle/10665/67462/WHO_SE_68.3.pdf>
"smallpox"

