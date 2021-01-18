#' Probability calculations for Weibull count models
#'
#' Probability computations for the univariate Weibull count process. Several
#' methods are provided.
#' \code{dWeibullCount} computes probabilities.
#'
#' Argument \code{method} can be used to specify the desired method, as follows:
#' \describe{
#' \item{\code{"series_mat"}}{- series expansion using matrix techniques,}
#' \item{\code{"series_acc"}}{- Euler-van Wijngaarden accelerated series expansion (default),}
#' \item{\code{"conv_direc"t}}{- direct convolution method of section 2,}
#' \item{\code{"conv_naive"}}{- naive convolurion described in section 3.1,}
#' \item{\code{"conv_dePril"}}{- dePril convolution described in section 3.2.}
#' }
#' The arguments have sensible default values.
#'
#' @param scale numeric (length 1), scale parameter of the Weibull count.
#' @param shape numeric (length 1), shape parameter of the Weibull count.
#' @param x integer (vector), the desired count values.
#' @param time double, length of the observation window (defaults to 1).
#' @param log logical, if TRUE, the log of the probability will be returned.
#' @param method character, one of the available methods, see details.
#' @param conv_steps numeric, number of steps used for the extrapolation.
#' @param conv_extrap logical, should Richardson extrappolation be applied ?
#' @param series_terms numeric, number of terms in the series expansion.
#' @param series_acc_niter numeric, number of iterations in the
#'     Euler-van Wijngaarden algorithm.
#' @param series_acc_eps numeric, tolerance of convergence in the
#'     Euler-van Wijngaarden algorithm.
#' @return for \code{dWeibullCount}, a vector of probabilities
#'     \eqn{P(x(i)), i = 1, \dots n}, where \eqn{n =} \code{length(x)}.
#' @export
dWeibullCount <- function(x, shape, scale,
                          method = c("series_acc", "series_mat",
                              "conv_direct", "conv_naive", "conv_dePril"),
                          time = 1, log = FALSE, conv_steps = 100,
                          conv_extrap = TRUE, series_terms = 50,
                          series_acc_niter = 300, series_acc_eps = 1e-10) {

    method <- match.arg(method)

     .dist <- function(i, scale, shape)
        list(scale = scale[i], shape = shape[i])

    vec <- FALSE
    if (length(x) == 1)
        distPars <- list(scale = scale, shape = shape)
    else {
        if (length(shape) == 1 & length(scale) == 1)
            distPars <- list(scale = scale, shape = shape)
        else if (length(shape) == 1 & length(scale) > 1) {
            shape <- rep(shape, length(scale))
            distPars <- lapply(1:length(shape), .dist,
                               scale = scale, shape = shape)
            vec <- TRUE
        } else if (length(shape) > 1 & length(scale) == 1) {
            scale <- rep(scale, length(shape))
            distPars <- lapply(1:length(shape), .dist,
                               scale = scale, shape = shape)
            vec <- TRUE
        } else {
            distPars <- lapply(1:length(shape), .dist,
                               scale = scale, shape = shape)
            vec <- TRUE
        }
    }
    dist <- "weibull"

    switch(method,
           "conv_direct" = {
              if (vec)
                   return(dCount_allProbs_vec_bi(x, distPars, dist,
                                                 conv_steps, time, conv_extrap,
                                                 log)
                          )
               else
                   return(dCount_allProbs_bi(x, distPars, dist,
                                             conv_steps, time, conv_extrap, log)
                          )
           },
           "conv_naive" = {
               if (vec)
                   return(dCount_naive_vec_bi(x, distPars, dist,
                                              conv_steps, time, conv_extrap,
                                              log))
               else
                   return(dCount_naive_bi(x, distPars, dist,
                                          conv_steps, time, conv_extrap, log))
           },
           "conv_dePril" = {
               if (vec)
                   return(dCount_dePril_vec_bi(x, distPars, dist,
                                               conv_steps, time, conv_extrap,
                                               log))
               else
                   return(dCount_dePril_bi(x, distPars, dist,
                                           conv_steps, time, conv_extrap,
                                           log))
           },
           "series_mat" = {
               if (vec)
                   return(dWeibullCount_mat_vec(x, shape, scale, time, log,
                                                series_terms))
               else
                   return(dWeibullCount_mat(x, shape, scale, time, log,
                                            series_terms))
           },
           "series_acc" = {
               if (vec)
                   return(dWeibullCount_acc_vec(x, shape, scale, time, log,
                                                series_terms, series_acc_niter,
                                                series_acc_eps)
                          )
               else
                   return(dWeibullCount_acc(x, shape, scale, time, log,
                                            series_terms, series_acc_niter,
                                            series_acc_eps)
                          )
           }
           )
}



#' % Log-likelihood of the the Weibull count process
#' % Compute the log-likelihood of the univariate Weibull count process. Several
#' % methods are provided.
#'
#' \code{dWeibullCount_loglik} computes the log-likelihood.
#'
#' @param na.rm logical, if TRUE \code{NA}'s (produced by taking the log of very
#'     small probabilities) will be replaced by the smallest allowed probaility;
#'     default is \code{TRUE}.
#' @param weights numeric, vector of weights to apply. If \code{NULL}, a vector
#'     of one's will be applied.
#' @inheritParams dWeibullCount
#'
#' @return for \code{dWeibullCount_loglik},
#'     a double, the log-likelihood of the count process.
#' @rdname dWeibullCount
#' @export
dWeibullCount_loglik <- function(x, shape, scale,
                                 method = c("series_acc", "series_mat",
                                     "conv_direct", "conv_naive", "conv_dePril"),
                                 time = 1, na.rm = TRUE, conv_steps = 100,
                                 conv_extrap = TRUE, series_terms = 50,
                                 series_acc_niter = 300,
                                 series_acc_eps = 1e-10, weights = NULL) {
    if (is.null(weights))
        weights <- rep(1, length(x))

    method <- match.arg(method)

   pbs <- dWeibullCount(x, shape, scale, method,
                         time = time, log = TRUE, conv_steps = conv_steps,
                         conv_extrap = conv_extrap, series_terms = series_terms,
                         series_acc_niter = series_acc_niter,
                         series_acc_eps = series_acc_eps)

    pbs <- pbs * weights
     if (na.rm)
         pbs[is.na(pbs)] <- .logNaReplace()

     sum(pbs)
}


#' % Weibull count expected value and variance
#' % Expected value and variance of the weibull count process
#'
#' \code{evWeibullCount} computes the expected value and variance.
#'
#' @param xmax unsigned integer, maximum count to be used.
#' @inheritParams dWeibullCount
#' @return for \code{evWeibullCount}, a list with components:
#'     \item{ExpectedValue}{expected value,}
#'     \item{Variance}{variance.}
#' @rdname dWeibullCount
#' @export
evWeibullCount <- function(xmax, shape, scale,
                           method = c("series_acc", "series_mat",
                               "conv_direct", "conv_naive", "conv_dePril"),
                           time = 1, conv_steps = 100,
                           conv_extrap = TRUE, series_terms = 50,
                           series_acc_niter = 300,
                           series_acc_eps = 1e-10) {

    method <- match.arg(method)

    x <- 1:xmax
    px <- dWeibullCount(x, shape, scale, method, time, FALSE,
                        conv_steps, conv_extrap, series_terms,
                        series_acc_niter, series_acc_eps)
    
    ## adjust probability 
    px[px < 0] <- 0
    px <- zapsmall(px)
    ##px <- px / sum(px)
    
    ev <- sum(x * px)
    ev2 <- sum(x^2 * px)
    var <- ev2 - ev^2
    if (var < 0) {
        var <- NA
        warning("failed to compute variance! NA returned")
    }
    
    list(ExpectedValue = ev, Variance = var)
}
