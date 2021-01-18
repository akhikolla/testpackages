#' Numerically stable log of the inv_logit function
#' @keywords internal
log_inv_logit <- function(u) {
    ifelse(u < 0, u - log1p(exp(u)), - log1p(exp(-u)))
}
