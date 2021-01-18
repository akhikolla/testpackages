#' Probability calculations for Weibull-gamma count models
#'
#' Probability computations for the univariate Weibull-gamma count processes. Several
#' methods are provided.
#' \code{dWeibullgammaCount} computes probabilities.
#'
#' The desired method can be specified by argument \code{method}, as follows:
#' \describe{
#' \item{\code{"series_mat"}}{series expansion using matrix techniques.}
#' \item{\code{"series_acc"}}{Euler-van Wijngaarden accelerated series expansion.}
#' }
#' The arguments have sensible default values.
#'
#' @param shapeGam,scaleGam numeric shape and scale parameters of the gamma
#'     heterogeity function.
#' @param shape numeric (length 1), shape parameter of the Weibull count.
#' @param x integer (vector), the desired count values.
#' @param Xcovar matrix the regressor values. Should have the same number of
#'     rows as \code{length(x)}. If NULL, no regression will be considered.
#' @param beta numeric regression coefficients. If NULL, no regression will be
#'     considered.
#' @param time double, length of the observation window (defaults to 1).
#' @param log logical, if TRUE, the log of the probability will be returned.
#' @param method character one of the available methods. See details.
#' @param series_terms numeric number of terms in the series expansion.
#' @param series_acc_niter numeric number of iteration in the Euler-van
#'     Wijngaarden algorithm.
#' @param series_acc_eps numeric tolerance of convergence in the Euler-van
#'     Wijngaarden algorithm.
#'
#' @return
#'     for \code{dWeibullgammaCount}, a vector of probabilities
#'     \eqn{P(x(i)), i = 1, \dots n} where \code{n = length(x)}.
#' @export
dWeibullgammaCount <- function(x, shape, shapeGam, scaleGam,
                               Xcovar = NULL, beta = NULL,
                               method = c("series_acc", "series_mat"),
                               time = 1, log = FALSE,  series_terms = 50,
                               series_acc_niter = 300, series_acc_eps = 1e-10) {

    method <- match.arg(method)
    ## -------------- rescale parameters if needed
    if (is.null(Xcovar) & is.null(beta))
        noReg <- TRUE
    else if (!is.null(Xcovar) & !is.null(beta)) {
        noReg <- FALSE

        if (length(x) != nrow(Xcovar))
            stop("length x should be the same as nrow Xcovar !")

        if (ncol(Xcovar) != length(beta))
            stop("check dimension Xcovar and beta !")
    } else
        stop("check inputs Xcovar and/or beta !")

    switch(method,
           "series_mat" = {
               if (noReg) {
                   if (length(x) == length(shape))
                       return(
                           as.numeric(
                               dWeibullgammaCount_mat_vec(x, shape,
                                                          shapeGam, scaleGam,
                                                          time, log, series_terms)
                           )
                       )
                   else
                       return(
                           as.numeric(
                               dWeibullgammaCount_mat(x, shape, shapeGam, scaleGam,
                                                      time, log, series_terms)
                           )
                       )
               } else {
                   if (length(x) == length(shape))
                       as.numeric(
                           dWeibullgammaCount_mat_Covariates_vec(x, shape,
                                                                 shapeGam,
                                                                 scaleGam,
                                                                 Xcovar, beta,
                                                                 time, log,
                                                                 series_terms)
                       )
                  else
                       return(
                           as.numeric(
                               dWeibullgammaCount_mat_Covariates(x, shape, shapeGam,
                                                                 scaleGam,
                                                                 Xcovar, beta,
                                                                 time, log,
                                                                 series_terms)
                           )
                       )
               }
           },
           "series_acc" = {
               if (noReg) {
                   if (length(x) == length(shape))
                       return(
                           as.numeric(
                               dWeibullgammaCount_acc_vec(x, shape,
                                                          shapeGam, scaleGam,
                                                          time, log, series_terms,
                                                          series_acc_niter,
                                                          series_acc_eps)
                           )
                       )
                   else
                       return(
                           as.numeric(
                               dWeibullgammaCount_acc(x, shape, shapeGam,
                                                      scaleGam,
                                                      time, log, series_terms,
                                                      series_acc_niter,
                                                      series_acc_eps)
                           )
                       )
               } else {
                   if (length(x) == length(shape))
                       return(
                           as.numeric(
                               dWeibullgammaCount_acc_Covariates_vec(
                                   x, shape,
                                   shapeGam,
                                   scaleGam,
                                   Xcovar, beta,
                                   time, log,
                                   series_terms,
                                   series_acc_niter,
                                   series_acc_eps)
                           )
                       )
                   else
                       return(
                           as.numeric(
                               dWeibullgammaCount_acc_Covariates(
                                   x, shape, shapeGam,
                                   scaleGam, Xcovar, beta,
                                   time, log,
                                   series_terms,
                                   series_acc_niter,
                                   series_acc_eps)
                           )
                       )
               }
           }
           )
}



#' % wrapper to the weibull-gamma count log-likelihood
#'
#' \code{dWeibullgammaCount_loglik} computes the log-likelihood.
#'
#' @param na.rm logical, if TRUE, \code{NA}'s (produced by taking the log of
#'     very small probabilities) will be replaced by the smallest allowed
#'     probaility; default = \code{TRUE}.
#' @param weights numeric, vector of weights to apply. If \code{NULL}, one will
#'     be applied.
#' @inheritParams dWeibullgammaCount
#'
#' @return
#'     for \code{dWeibullgammaCount_loglik}, double,
#'     log-likelihood of the count process
#' @rdname dWeibullgammaCount
#' @export
dWeibullgammaCount_loglik <- function(x, shape, shapeGam, scaleGam,
                                      Xcovar = NULL, beta = NULL,
                                      method = c("series_acc", "series_mat"),
                                      time = 1, na.rm = TRUE, series_terms = 50,
                                      series_acc_niter = 300,
                                      series_acc_eps = 1e-10, weights = NULL) {
    if (is.null(weights))
        weights <- rep(1, length(x))

    method <- match.arg(method)

       pbs <- dWeibullgammaCount(x, shape, shapeGam, scaleGam,
                              Xcovar, beta, method,
                              time = time, log = TRUE,
                              series_terms = series_terms,
                              series_acc_niter = series_acc_niter,
                              series_acc_eps = series_acc_eps)

     pbs <- as.numeric(pbs) * weights
     if (na.rm)
         pbs[is.na(pbs)] <- .logNaReplace()

     sum(pbs)
}


#' Weibull-gamma count expected value and variance
#'
#' \code{evWeibullgammaCount} computes the expected value and variance.
#'
#' @param xmax unsigned integer, maximum count to be used.
#' @inheritParams dWeibullgammaCount
#' @return
#'     for \code{evWeibullgammaCount}, a list with components
#'     \code{"ExpectedValue"} and \code{"Variance"}.
#' @rdname dWeibullgammaCount
#' @export
evWeibullgammaCount <- function(xmax, shape, shapeGam, scaleGam,
                                Xcovar = NULL, beta = NULL,
                                method = c("series_acc", "series_mat"),
                                time = 1, series_terms = 50,
                                series_acc_niter = 300, series_acc_eps = 1e-10) {

    method <- match.arg(method)
    ## make sure we have a one row object corresponding to individual i
    if (!is.null(Xcovar)) {
        stopifnot(nrow(Xcovar) == 1)
        Xcovar <- matrix(Xcovar, nrow = xmax, ncol = ncol(Xcovar),
                         byrow = TRUE)
    }

    x <- 1:xmax
    px <- dWeibullgammaCount(x, shape, shapeGam, scaleGam, Xcovar, beta,
                             method, time, FALSE,  series_terms,
                             series_acc_niter, series_acc_eps)

    ev <- sum(x * px)
    ev2 <- sum(x^2 * px)
    var <- ev2 - ev^2
    list(ExpectedValue = ev, Variance = var)
}
