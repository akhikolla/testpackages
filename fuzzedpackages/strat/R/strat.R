#' Stratification index.
#'
#' \code{strat} computes the stratification index proposed in Zhou (2012). When
#' \code{group} is specified, it also returns between-group and within-group
#' components of the overall stratification.
#'
#' @param outcome A numeric vector of outcome.
#' @param strata A vector of \code{length(outcome)} indicating strata
#'   membership. The elements are coerced to factors by
#'   \code{\link[base]{factor}}.
#' @param weights An optional vector of weights.
#' @param ordered Logical. If \code{TRUE} strata are pre-ordered ascendingly.
#' @param group An optional grouping factor. If specified, \code{strat} also
#'   returns between-group and within-group components of the overall
#'   stratification.
#' @return An object of class \code{strat}. \item{overall}{a vector of two,
#'   giving computed stratification index and approximate standard error.}
#'   \item{strata_info}{a data frame of stratum-specific information, including
#'   name, population share, and average percentile
#'   rank.}\item{decomposition}{between-group and within-group components of the
#'   overall stratification.}\item{within_group}{within-group indices of
#'   stratification by group.}
#' @export
#' @references Zhou, Xiang. 2012. "A Nonparametric Index of Stratification."
#'   Sociological Methodology, 42(1): 365-389.
#' @examples
#' s <- with(cpsmarch2015, strat(income, big_class,
#'  weights = weight, group = education))
#' print(s, digits = 4)
#' print(s$strata_info, digits = 4)
#' print(s$within_group, digits = 4)
strat <- function(outcome, strata, weights = NULL,
    ordered = FALSE, group = NULL) {

    if (!is.logical(ordered) || is.na(ordered))
      stop("ordered has to be a valid logical scalar")

    if (anyNA(group))
      stop("group contains missing values")

    group_name <- deparse(substitute(group))

    input <- srank(outcome, strata, weights, group)
    input_df <- input$raw
    strata_info <- input$summary
    input_df <- merge(input_df, strata_info, by.x = "strata",
        by.y = "strata")
    input_df$sort_by <- ordered * as.numeric(input_df$strata) +
        (1 - ordered) * input_df$s_prank
    input_df <- input_df[order(input_df$sort_by), ]

    if (length(unique(input_df$group)) == 1) {
        s <- strat_cpp(input_df$prank, input_df$sort_by,
            input_df$weights)
        decomp <- NULL
        within <- NULL
    } else {
        group_num <- as.numeric(input_df$group) - 1
        s <- strat_cpp_by(input_df$prank, input_df$sort_by,
            input_df$weights, group_num)
        decomp <- matrix(c(s$weight_w, s$weight_b,
            s$strat_w, s$strat_b), 2, 2, dimnames = list(c(paste("within",
            group_name), paste("between", group_name)),
            c("weight", "strat")))
        within <- data.frame(group = levels(input_df$group),
            weight = s$weight_by, strat = s$strat_by)
        names(within)[1] <- group_name
    }

    # approximate standard error (Goodman & Kruskal 1963)
    arg <- s$deno/(1 - (s$strat^2))/nrow(input_df)
    std_error <- 1/sqrt(arg)

    overall <- c(strat = s$strat, std_error = std_error)

    # return a list
    out <- list(overall = overall, strata_info = strata_info,
        decomposition = decomp, within_group = within)
    class(out) <- c("strat", "list")
    out
}




