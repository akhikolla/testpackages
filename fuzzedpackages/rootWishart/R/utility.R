# Error function
#' @importFrom stats pnorm
erf <- function(x) 2*pnorm(sqrt(2)*x) - 1

# Special cases for single Wishart----
F22 <- function(x) {
    1 - sqrt(x*pi*0.5)*exp(-0.5*x)*erf(sqrt(0.5*x)) - exp(-x)
}
F25 <- function(x) {
    1 - exp(-x)*(2*exp(0.5*x)*x^2 + x^2 + 6*x + 6)/6
}
F33 <- function(x) {
    exp(-1.5*x)*(exp(0.5*x) * (exp(x) - x - 1) * erf(sqrt(0.5*x)) -
                     sqrt(2*x/pi)*(exp(x)*(x - 1) + 1))
}
F44 <- function(x) {
    exp(-2*x)*(sqrt(2)*(4*exp(2*x) -
                            exp(x)*(x^3 + 2*x^2 + 2*x + 8) + 2*x + 4) -
                   sqrt(pi*x)*exp(0.5*x)*(exp(x)*(x^2 - 4*x + 6) -
                                              2*x - 6)*erf(sqrt(0.5*x)))/sqrt(32)
}

# Check if we have an integer----
isWhole <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# Adaptive selection of precision type----
precSingleWishart_bool <- function(n_min, n_max) {
    multiPrec <- FALSE

    if (n_max > 17) multiPrec <- TRUE else {
        if (n_min >= 15) multiPrec <- TRUE
        if (n_min == 14 & (n_max %in% 15:16)) multiPrec <- TRUE

    }
    return(multiPrec)
}

precDoubleWishart_bool <- function(p, n, m) {
    multiPrec <- FALSE

    if (p > 12 || m > 250 || n > 40) multiPrec <- TRUE
    if (abs(p - n) > 2) multiPrec <- TRUE

    return(multiPrec)
}
