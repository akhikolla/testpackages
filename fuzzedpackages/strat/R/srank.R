#' Ranking strata.
#'
#' Ranking strata according to the average percentile rank of members in each
#' stratum.
#'
#' @inheritParams strat
#' @param group An optional grouping factor.
#' @return An object of class \code{srank}. \item{raw}{a data frame
#'   consisting of complete cases of all inputs.} \item{summary}{a data frame
#'   of stratum-specific information, including name, population
#'   share, and average percentile rank.}
#' @export
#' @examples
#' strata_info <- with(cpsmarch2015, srank(income, big_class,
#'  weights = weight, group = education))
#' print(strata_info, digits = 3)
srank <- function(outcome, strata, weights = NULL,
    group = NULL) {

    # data input
    input <- clean(outcome, strata, weights = weights, group = group)

    # ranking strata
    tmp_fun <- function(x) c(share = sum(x$weights)/sum(input$weights),
        s_prank = stats::weighted.mean(x$prank, x$weights))
    out <- t(vapply(split(input, input$strata),
        tmp_fun, numeric(2)))

    # output
    strata_info <- data.frame(strata = rownames(out), out)
    rownames(strata_info) <- NULL

    out <- list(raw = input, summary = strata_info)
    class(out) <- c("srank", "list")
    out
}
