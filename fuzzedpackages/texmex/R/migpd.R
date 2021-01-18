#' Fit multiple independent generalized Pareto models
#'
#' Fit multiple independent generalized Pareto models as the first step of
#' conditional multivariate extreme values modelling following the approach of
#' Heffernan and Tawn, 2004.
#'
#' The parameters in the generalized Pareto distribution are estimated for each
#' column of the data in turn, independently of all other columns. Note,
#' covariate modelling of GPD parameters is not supported.
#'
#' Maximum likelihood estimation often fails with generalized Pareto
#' distributions because of the likelihood becoming flat (see, for example,
#' Hosking et al, 1985).  Therefore the function allows penalized likelihood
#' estimation, which is the same as maximum a posteriori estimation from a
#' Bayesian point of view.
#'
#' By default quadratic penalization is used, corresponding to using a Gaussian
#' prior. If no genuine prior information is available, the following argument
#' can be used. If xi = -1, the generalized Pareto distribution corresponds to
#' the uniform distribution, and if xi is 1 or greater, the expectation is
#' infinite. Thefore, xi is likely to fall in the region (-1, 1). A Gaussian
#' distribution centred at zero and with standard deviation 0.5 will have
#' little mass outside of (-1, 1) and so will often be a reasonable prior for
#' xi. For log(sigma) a Gaussian distribution, centred at zero and with
#' standard deviation 100 will often be vague.  If a Gaussian penalty is
#' specified but no parameters are given, the function will assume such
#' indpendent priors.
#'
#' Note that internally the function works with log(sigma), not sigma. The
#' reasons are that quadratic penalization makes more sense for phi=log(sigma)
#' than for sigma (because the distribution of log(sigma) will be more nearly
#' symmetric), and because it was found to stabilize computations.
#'
#' The associated \code{coef}, \code{print} and \code{summary} functions
#' exponentiate the log(sigma) parameter to return results on the expected
#' scale. If you are accessesing the parameters directly, however, take care to
#' be sure what scale the results are on.
#'
#' Threshold selection can be carried out with the help of functions
#' \code{\link{mrl}} and \code{\link{gpdRangeFit}}.
#'
#' @aliases migpd plot.migpd ggplot.migpd
#' @note You are encourage to use the \code{mqu} argument and not \code{mth}.
#'   If you use \code{mth}, the quantiles then need to be estimated. There
#'   are, at the time of writing, 9 methods of estimating quantiles build into
#'   the \code{quantile} function. Tiny differences can cause problems in
#'   later stages of the analysis if functions try to simulate in an area
#'   that is legitimate according to the numerical value of the threshold, but
#'   not according to the estimated quantile.
#' @param data A matrix or data.frame, each column of which is to be modelled.
#' @param mth Marginal thresholds. Thresholds above which to fit the models.
#' Only one of \code{mth} and \code{mqu} should be supplied. Length one (in
#' which case a common threshold is used) or length equal to the number of
#' columns of \code{data} (in which case values correspond to thresholds for
#' each of the columns respectively).
#' @param mqu Marginal quantiles. Quantiles above which to fit the models. Only
#' one of \code{mth} and \code{mqu} should be supplied. Length as for
#' \code{mth} above.
#' @param penalty How the likelihood should be penalized. Defaults to
#' "gaussian". See documentation for \code{\link{evm}}.
#' @param maxit The maximum number of iterations to be used by the optimizer.
#' @param trace Whether or not to tell the user how the optimizer is getting
#' on. The argument is passed into \code{\link{optim}} -- see the help for that
#' function.
#' @param verbose Controls whether or not the function prints to screen every
#' time it fits a model. Defaults to FALSE.
#' @param priorParameters Only used if \code{penalty = 'gaussian'}. A named
#' list, each element of which contains two components: the first should be a
#' vector of length 2 corresponding to the location of the Gaussian
#' distribution; the second should be 2x2 matrix corresponding to the
#' covariance matrix of the distribution. The names should match the names of
#' the columns of \code{data}. If not provided, it defaults to independent
#' priors being centred at zero, with variance 10000 for log(sigma) and 0.25
#' for xi. See the details section.
#' @param cov String, passed through to \code{evm}: how to estimate the covariance.
#'   Defaults to \code{cov = "observed"}.
#' @param family An object of class "texmexFamily". Should be either
#'   \code{family = gpd} or \code{family = cgpd} and defaults to the first of those.
#' @param x Object of class \code{migpd} as returned by function \code{migpd}.
#' @param main Character vector of length four: titles for plots produced by
#' \code{plot} and \code{ggplot} methods.
#' @param xlab As \code{main} but for x-axes labels.
#' @param nsim Number of simulations on which to base tolerance envelopes in
#' \code{plot} and \code{ggplot} methods.
#' @param alpha Significance level for tolerance and confidence intervals in
#' \code{plot} and \code{ggplot} methods.
#' @param mapping,environment Further arguments to ggplot method.
#' @param ... Further arguments to be passed to methods.
#' @return An object of class "migpd". There are \code{coef}, \code{print},
#' \code{plot}, \code{ggplot} and \code{summary} functions available.
#' @author Harry Southworth
#' @seealso \code{\link{mex}}, \code{\link{mexDependence}},
#' \code{\link{bootmex}}, \code{\link{predict.mex}}, \code{\link{gpdRangeFit}},
#' \code{\link{mrl}}
#' @references J. E. Heffernan and J. A. Tawn, A conditional approach for
#' multivariate extreme values, Journal of the Royal Statistical society B, 66,
#' 497 -- 546, 2004
#'
#' J. R. M. Hosking and J. R. Wallis, Parameter and quantile estimation for the
#' genralized Pareto distribution, Technometrics, 29, 339 -- 349, 1987
#' @keywords models multivariate
#' @examples
#'
#' \donttest{
#' mygpd <- migpd(winter, mqu=.7, penalty = "none")
#' mygpd
#' summary(mygpd)
#' plot(mygpd)
#' g <- ggplot(mygpd)
#' }
#'
#' @export migpd
migpd <-
function (data, mth, mqu, penalty = "gaussian", maxit = 10000,
   trace = 0, verbose=FALSE, priorParameters = NULL, cov = "observed",
   family = gpd){
   theCall <- match.call()
   if(is.null(colnames(data))){
    colnames(data) <- paste(rep("Column",ncol(data)),1:ncol(data),sep="")
  }
   d <- dim(data)[2]
   if (missing(mth) & missing(mqu))
       stop("you must provide one of mth or mqu")
   if (!missing(mth) & !missing(mqu))
       stop("you must provide precisely one of mth or mqu")
   if (!(family$name %in% c("GPD", "CGPD"))){
     stop("family should be either gpd or cgpd")
   }
   if (!missing(mth))
       mth <- rep(mth, length = d)
   if (!missing(mqu))
       mqu <- rep(mqu, length = d)
   if (missing(mqu))
       mqu <- sapply(1:d, function(i, x, mth) 1 - mean(x[, i] >
           mth[i]), x = data, mth = mth)
   if (missing(mth))
       mth <- sapply(1:d, function(i, x, prob) quantile(x[, i],
           prob = prob[i]), x = data, prob = mqu)
   if (penalty %in% c("quadratic", "gaussian") & is.null(priorParameters)) {
       gp = list(c(0, 0), matrix(c(100^2, 0, 0, 0.25), nrow = 2))
       priorParameters <- vector("list", length = length(mth))
       for (i in 1:length(mth)) priorParameters[[i]] <- gp
       names(priorParameters) <- dimnames(data)[[2]]
   } else if (penalty %in% c("quadratic", "gaussian")) {
       nm <- names(priorParameters)
       if (is.null(nm)){
           stop("priorParameters must be a named list")
       } else if (any(!is.element(nm, dimnames(data)[[2]]))){
           stop("the names of priorParameters must match the column names of the data")
       }
   }

   wrapgpd <-function(i, x, mth, penalty, maxit, verbose, trace, priorParameters) {
       if (verbose)
           cat("Fitting model", i, "\n")
       if (!is.null(priorParameters))
           priorParameters <- priorParameters[[(1:length(priorParameters))[names(priorParameters) ==
               dimnames(x)[[2]][i]]]]
       x <- c(x[, i])
       mth <- mth[i]

       evm(x, th=mth, penalty=penalty, priorParameters=priorParameters,
           maxit=maxit, trace=trace, cov = cov, family = family)
       }

   modlist <- lapply(1:d, wrapgpd, x=data, penalty=penalty, mth=mth, verbose=verbose,
                        priorParameters=priorParameters,maxit=maxit,trace=trace)
   if (length(dimnames(data)[[2]]) == dim(data)[[2]]){
       names(modlist) <- dimnames(data)[[2]]
   }
   names(mth) <- names(mqu) <- dimnames(data)[[2]]
   res <- list(call = theCall, models = modlist, data = data,
       mth = mth, mqu = mqu, penalty = penalty, priorParameters = priorParameters)
   oldClass(res) <- "migpd"
   invisible(res)
}

