## Utility functions used in the rocTree package


#' Used to calculate the empirical density function from ecdf
#'
#' @noRd
#' @keywords internal
fecdf <- function(x) rank(x, ties.method = "max", na.last = "keep") / length(x)


#' Smoothing kernel
#'
#' @keywords internal
#' @noRd
K1 <- function(u) {
  0.75 * (1 - u ^ 2) * (abs(u) < 1)
}

#' Smoothing kernel
#' This doesn't adjust for boundary condition and can cause errors
#' @keywords internal
#' @noRd
K2 <- function(s, vec, h) {
    if (is.na(s)) return(rep(NA, length(vec)))
    if (s < h) return(K1((h - vec) / h))
    ## if (s < h) return(Kq((s - vec) / h, s / h))
    else return(K1((s - vec) / h))
}

#' Smoothing kernel
#' @keywords internal
#' @noRd
Kq <- function(x, q) {
    sigk1 <- sqrt(0.2)
    2 / (q + 1) * K1(2 / (q + 1) * (x - (q - 1) / 2)) *
        (1 + ((q - 1) / (q + 1) / sigk1) ^ 2 + 2 / sigk1 ^ 2 * (1 - q) / (1 + q) ^ 2 * x)
}

#' Smoothing kernels for hazard estimation
#' @keywords internal
#' @noRd
K3 <- function(s, vec, h) {
    if (is.na(s)) return(rep(NA, length(vec)))
    if (s < h) return(Kq((s - vec) / h, s / h))
    else return(K1((s - vec) / h))
}


#' Check if input is a wholenumber
#' @keywords internal
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
