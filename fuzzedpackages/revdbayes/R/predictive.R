# =========================== predict.evpost ===========================

#' Predictive inference for the largest value observed in N years.
#'
#' \code{predict} method for class "evpost".  Performs predictive inference
#' about the largest value to be observed over a future time period of
#' N years.  Predictive inferences accounts for uncertainty in model
#' parameters and for uncertainty owing to the variability of future
#' observations.
#'
#' @param object An object of class "evpost", a result of a call to
#'   \code{\link{rpost}} or \code{\link{rpost_rcpp}} with \code{model = "gev"},
#'   \code{model = "os"}, \code{model = "pp"} or \code{model == "bingp"}.
#'   Calling these functions after a call to \code{rpost} or \code{rpost_rcpp}
#'   with \code{model == "gp"} will produce an error, because inferences about
#'   the probability of threshold exceedance are required, in addition to the
#'   distribution of threshold excesses. The model is stored in
#'   \code{object$model}.
#' @param type A character vector.  Indicates which type of inference is
#'   required:
#' \itemize{
#'   \item "i" for predictive intervals,
#'   \item "p" for the predictive distribution function,
#'   \item "d" for the predictive density function,
#'   \item "q" for the predictive quantile function,
#'   \item "r" for random generation from the predictive distribution.
#' }
#' @param x A numeric vector or a matrix with \code{n_years} columns.
#'   The meaning of \code{x} depends on \code{type}.
#'   \itemize{
#'     \item{\code{type = "p"} or \code{type = "d"}:} \code{x} contains
#'       quantiles at which to evaluate the distribution or density function.
#'
#'       If \code{object$model == "bingp"} then no element of \code{x} can be
#'       less than the threshold \code{object$thresh}.
#'
#'       If \code{x} is not supplied then \code{n_year}-specific defaults are
#'       set: vectors of length \code{x_num} from the 0.1\% quantile to the
#'       99\% quantile, subject all values being greater than the threshold.
#'
#'     \item{\code{type = "q"}:} \code{x} contains probabilities in (0,1)
#'       at which to evaluate the quantile function.  Any values outside
#'       (0, 1) will be removed without warning.
#'
#'       If \code{object$model == "bingp"} then no element of \code{p} can
#'       correspond to a predictive quantile that is below the threshold,
#'       \code{object$thresh}.  That is, no element of \code{p} can be less
#'       than the value of \code{predict.evpost(object,}
#'       \code{type = "q", x = object$thresh)}.
#'
#'       If \code{x} is not supplied then a default value of
#'       \code{c(0.025, 0.25, 0.5, 0.75, 0.975)} is used.
#'     \item{\code{type = "i"} or \code{type = "r"}:} \code{x} is not relevant.
#'   }
#' @param x_num A numeric scalar.  If \code{type = "p"} or \code{type = "d"}
#'   and \code{x} is not supplied then \code{x_num} gives the number of values
#'   in \code{x} for each value in \code{n_years}.
#' @param n_years A numeric vector. Values of N.
#' @param npy A numeric scalar. The mean number of observations per year
#'   of data, after excluding any missing values, i.e. the number of
#'   non-missing observations divided by total number of years of non-missing
#'   data.
#'
#' If \code{rpost} or \code{rpost_rcpp} was called with
#' \code{model == "bingp"} then \code{npy} must either have been supplied
#' in that call or be supplied here.
#'
#' Otherwise, a default value will be assumed if \code{npy} is not supplied,
#' based on the value of \code{model} in the call to \code{rpost} or
#' \code{rpost_rcpp}:
#' \itemize{
#'   \item{\code{model = "gev"}:} \code{npy} = 1, i.e. the data were
#'     annual maxima so the block size is one year.
#'   \item{\code{model = "os"}:} \code{npy} = 1, i.e. the data were
#'     annual order statistics so the block size is one year.
#'   \item{\code{model = "pp"}:}
#'     \code{npy} = \code{length(x$data)} / \code{object$noy},
#'     i.e. the value of \code{noy} used in the call to \code{\link{rpost}}
#'     or \code{\link{rpost_rcpp}} is equated to a block size of one year.
#' }
#' If \code{npy} is supplied twice then the value supplied here will be
#' used and a warning given.
#' @param level A numeric vector of values in (0, 100).
#'   Only relevant when \code{type = "i"}.
#'   Levels of predictive intervals for the largest value observed in
#'   N years, i.e. level\% predictive intervals are returned.
#' @param hpd A logical scalar.
#'   Only relevant when \code{type = "i"}.
#'
#'   If \code{hpd = FALSE} then the interval is
#'   equi-tailed, with its limits produced by
#'   \code{predict.evpost(}\code{object, type ="q", x = p)},
#'   where \code{p = c((1-level/100)/2,} \code{(1+level/100)/2)}.
#'
#'   If \code{hpd = TRUE} then, in addition to the equi-tailed interval,
#'   the shortest possible level\% interval is calculated.
#'   If the predictive distribution is unimodal then this
#'   is a highest predictive density (HPD) interval.
#' @param lower_tail A logical scalar.
#'   Only relevant when \code{type = "p"} or \code{type = "q"}.
#'   If TRUE (default), (output or input) probabilities are
#'   \eqn{P[X \leq x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @param log A logical scalar.  Only relevant when \code{type = "d"}.
#'   If TRUE the log-density is returned.
#' @param big_q A numeric scalar.  Only relevant when \code{type = "q"}.
#'   An initial upper bound for the desired quantiles to be passed to
#'   \code{\link[stats]{uniroot}} (its argument \code{upper}) in the
#'   search for the predictive quantiles.  If this is not sufficiently large
#'   then it is increased until it does provide an upper bound.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details Inferences about future extreme observations are integrated over
#'   the posterior distribution of the model parameters, thereby accounting
#'   for uncertainty in model parameters and uncertainty owing to the
#'   variability of future observations.  In practice the integrals involved
#'   are estimated using an empirical mean over the posterior sample.
#'   See, for example,
#'   \href{https://doi.org/10.1007/978-1-4471-3675-0_9}{Coles (2001),
#'   chapter 9},
#'   \href{https://doi.org/10.1201/b19721}{Stephenson (2016)}
#'   or
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   for details.
#'   See also the vignette
#'   \href{https://CRAN.R-project.org/package=revdbayes}{Posterior Predictive Extreme Value Inference}
#'
#'   \strong{GEV / OS / PP}.
#'   If \code{model = "gev"}, \code{model = "os"} or \code{model = "pp"}
#'   in the call to \code{\link{rpost}} or \code{\link{rpost_rcpp}}
#'   we first calculate the number of blocks \eqn{b} in \code{n_years} years.
#'   To calculate the density function or distribution function of the maximum
#'   over \code{n_years} we call \code{\link{dgev}} or \code{\link{pgev}}
#'   with \code{m} = \eqn{b}.
#'
#'   \itemize{
#'     \item{\code{type = "p"}.} We calculate using \code{\link{pgev}}
#'     the GEV distribution function at \code{q} for each of the posterior
#'     samples of the location, scale and shape parameters.  Then we take
#'     the mean of these values.
#'
#'     \item{\code{type = "d"}.} We calculate using \code{\link{dgev}}
#'     the GEV density function at \code{x} for each of the posterior samples
#'     of the location, scale and shape parameters.  Then we take the
#'     mean of these values.
#'
#'     \item{\code{type = "q"}.} We solve numerically
#'     \code{predict.evpost(object, type = "p", x = q)} = \code{p[i]}
#'     numerically for \code{q} for each element \code{p[i]} of \code{p}.
#'
#'     \item{\code{type = "i"}.} If \code{hpd = FALSE} then the interval is
#'     equi-tailed, equal to \code{predict.evpost()} \code{object, type ="q", x = p)},
#'     where \code{p = c((1-level/100)/2,} \code{(1+level/100)/2)}.
#'     If \code{hpd = TRUE} then, in addition, we perform a
#'     numerical minimisation of the length of level\% intervals, after
#'     approximating the predictive quantile function using monotonic
#'     cubic splines, to reduce computing time.
#'
#'     \item{\code{type = "r"}.} For each simulated value of the GEV parameters
#'     at the \code{n_years} level of aggregation we simulate one value from
#'     this GEV distribution using \code{\link{rgev}}.  Thus, each sample
#'     from the predictive distribution is of a size equal to the size of
#'     the posterior sample.
#'   }
#'
#'   \strong{Binomial-GP}.  If \code{model = "bingp"} in the call to
#'   \code{\link{rpost}} or \code{\link{rpost_rcpp}} then we calculate the
#'   mean number of observations in \code{n_years} years, i.e.
#'   \code{npy * n_years}.
#'
#'   Following \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   let \eqn{M_N} be the largest value observed in \eqn{N} years,
#'   \eqn{m} = \code{npy * n_years} and \eqn{u} the threshold
#'   \code{object$thresh} used in the call to \code{rpost}
#'   or \code{rpost_rcpp}.
#'   For fixed values of \eqn{\theta = (p, \sigma, \xi)} the distribution
#'   function of \eqn{M_N} is given by \eqn{F(z, \theta)^m}, for
#'   \eqn{z \geq u}{z >= u}, where
#'   \deqn{F(z, \theta) = 1 - p [1 + \xi (x - u) / \sigma] ^ {-1/\xi}.}{%
#'         F(z, \theta) = 1 - p * [1 + \xi (x - u) / \sigma] ^ (-1/\xi).}
#'   The distribution function of \eqn{M_N} cannot be evaluated for
#'   \eqn{z < u} because no model has been supposed for observations below
#'   the threshold.
#'
#' \itemize{
#'   \item{\code{type = "p"}.} We calculate
#'     \eqn{F(z, \theta)^m} at \code{q} for each of the posterior samples
#'     \eqn{\theta}.  Then we take the mean of these values.
#'   \item{\code{type = "d"}.}  We calculate the density of of \eqn{M_n}, i.e.
#'     the derivative of \eqn{F(z, \theta)^m} with respect to \eqn{z} at
#'     \code{x} for each of the posterior samples \eqn{\theta}.  Then we take
#'     the mean of these values.
#'   \item{\code{type = "q"} and \code{type = "i"}} We perform calculations
#'     that are analogous to the GEV case above.  If \code{n_years} is very
#'     small and/or level is very close to 100 then a predictive interval
#'     may extend below the threshold.  In such cases \code{NA}s are returned
#'     (see \strong{Value} below).
#'   \item{\code{type = "r"}.}  For each simulated value of the bin-GP
#'     parameter we simulate from the distribution of \eqn{M_N} using the
#'     inversion method applied to the distribution function of \eqn{M_N} given
#'     above.  Occasionally a value below the threshold would need to be
#'     simulated.  If these instances a missing value code \code{NA} is
#'     returned. Thus, each sample from the predictive distribution is of a
#'     size equal to the size of the posterior sample, perhaps with a small
#'     number os \code{NA}s.
#'   }
#' @return An object of class "evpred", a list containing a subset of the
#'   following components:
#'     \item{type}{The argument \code{type} supplied to \code{predict.evpost}.
#'     Which of the following components are present depends \code{type}.}
#'     \item{x}{A matrix containing the argument \code{x} supplied to
#'       \code{predict.evpost}, or set within \code{predict.evpost} if \code{x}
#'       was not supplied, replicated to have \code{n_years} columns
#'       if necessary.
#'       Only present if \code{type} is \code{"p", "d"} or \code{"q"}.}
#'     \item{y}{The content of \code{y} depends on \code{type}:
#'     \itemize{
#'       \item{\code{type = "p", "d", "q"}:}  A matrix with the same
#'       dimensions as \code{x}.  Contains distribution function values
#'       (\code{type = "p"}), predictive density (\code{type = "d"})
#'       or quantiles (\code{type = "q"}).
#'       \item{\code{type = "r"}:} A numeric matrix with \code{length(n_years)}
#'       columns and number of rows equal to the size of the posterior sample.
#'       \item{\code{type = "i"}:} \code{y} is not present.
#'       }}
#'       \item{long}{A \code{length(n_years)*length(level)} by 4 numeric
#'         matrix containing the equi-tailed limits with columns:
#'         lower limit, upper limit, n_years, level.
#'         Only present if \code{type = "i"}.  If an interval extends below
#'         the threshold then \code{NA} is returned.}
#'       \item{short}{A matrix with the same structure as \code{long}
#'         containing the HPD limits.  Only present if \code{type = "i"}.
#'         Columns 1 and 2 contain \code{NA}s if \code{hpd = FALSE}
#'         or if the corresponding equi-tailed interval extends below
#'         the threshold.}
#'   The arguments \code{n_years, level, hpd, lower_tail, log} supplied
#'   to \code{predict.evpost} are also included, as is the argument \code{npy}
#'   supplied to, or set within, \code{predict.evpost} and
#'   the arguments \code{data} and \code{model} from the original call to
#'   \code{\link{rpost}} or \code{\link{rpost_rcpp}}.
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   Chapter 9: \url{https://doi.org/10.1007/978-1-4471-3675-0_9}
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
#' @references Stephenson, A. (2016). Bayesian Inference for Extreme Value
#'   Modelling. In \emph{Extreme Value Modeling and Risk Analysis: Methods and
#'   Applications}, edited by D. K. Dey and J. Yan, 257-80. London:
#'   Chapman and Hall. \url{https://doi.org/10.1201/b19721}
#' @seealso \code{\link{plot.evpred}} for the S3 \code{plot} method for
#'   objects of class \code{evpred}.
#' @seealso \code{\link{rpost}} or \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @examples
#' ### GEV
#' data(portpirie)
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' gevp  <- rpost_rcpp(n = 1000, model = "gev", prior = pn, data = portpirie)
#'
#' # Interval estimation
#' predict(gevp)$long
#' predict(gevp, hpd = TRUE)$short
#' # Density function
#' x <- 4:7
#' predict(gevp, type = "d", x = x)$y
#' plot(predict(gevp, type = "d", n_years = c(100, 1000)))
#' # Distribution function
#' predict(gevp, type = "p", x = x)$y
#' plot(predict(gevp, type = "p", n_years = c(100, 1000)))
#' # Quantiles
#' predict(gevp, type = "q", n_years = c(100, 1000))$y
#' # Random generation
#' plot(predict(gevp, type = "r"))
#'
#' ### Binomial-GP
#' data(gom)
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' bp <- set_bin_prior(prior = "jeffreys")
#' npy_gom <- length(gom)/105
#' bgpg <- rpost_rcpp(n = 1000, model = "bingp", prior = fp, thresh = u,
#'                    data = gom, bin_prior = bp)
#'
#' # Setting npy in call to predict.evpost()
#' predict(bgpg, npy = npy_gom)$long
#'
#' # Setting npy in call to rpost() or rpost_rcpp()
#' bgpg <- rpost_rcpp(n = 1000, model = "bingp", prior = fp, thresh = u,
#'                    data = gom, bin_prior = bp, npy = npy_gom)
#'
#' # Interval estimation
#' predict(bgpg)$long
#' predict(bgpg, hpd = TRUE)$short
#' # Density function
#' plot(predict(bgpg, type = "d", n_years = c(100, 1000)))
#' # Distribution function
#' plot(predict(bgpg, type = "p", n_years = c(100, 1000)))
#' # Quantiles
#' predict(bgpg, type = "q", n_years = c(100, 1000))$y
#' # Random generation
#' plot(predict(bgpg, type = "r"))
#' @export
predict.evpost <- function(object, type = c("i", "p", "d", "q", "r"), x = NULL,
                           x_num = 100, n_years = 100, npy = NULL, level = 95,
                           hpd = FALSE, lower_tail = TRUE, log = FALSE,
                           big_q = 1000, ...) {
  type <- match.arg(type)
  if (!inherits(object, "evpost")) {
    stop("object must be an evpost object produced by rpost() or rpost_rcpp()")
  }
  if (object$model == "gp") {
    stop("The model cannot be gp.  Use bingp instead.")
  }
  if (!(object$model %in% c("gev", "os", "pp", "bingp"))) {
    stop(paste("Predictive functions are not available for model = ''",
                object$model, "''", sep=""))
  }
  # Set the value of npy.
  npy <- set_npy(object = object, npy = npy)
  #
  n_y <- length(n_years)
  if (type == "i") {
    n_l <- length(level)
    ret_obj <- list()
    ret_obj$long <- ret_obj$short <- matrix(NA, ncol = 4, nrow = n_y * n_l)
    cnames <- c("lower", "upper", "n_years", "level")
    colnames(ret_obj$long) <- colnames(ret_obj$short) <- cnames
    if ((length(n_years) == 1 & length(level) == 1)) {
      temp <- ipred(object, n_years = n_years, npy = npy, level = level,
                    hpd = hpd, big_q = big_q)
      ret_obj$long[1, ] <- c(temp$long, n_years, level)
      ret_obj$short[1, ] <- c(temp$short, n_years, level)
    } else {
      k <- 1
      for (i in 1:n_y) {
        for (j in 1:n_l) {
        temp <- ipred(object, n_years = n_years[i], npy = npy,
                      level = level[j], hpd = hpd, big_q = big_q)
        ret_obj$long[k, ] <- c(temp$long, n_years[i], level[j])
        ret_obj$short[k, ] <- c(temp$short, n_years[i], level[j])
        k <- k + 1
        }
      }
    }
  }
  if (type == "q") {
    if (is.null(x)) {
      x <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    } else {
      x <- x[x > 0 & x < 1]
      if (length(x) == 0) {
        stop("For type = ``q'' values in x must be in (0,1)")
      }
    }
  }
  if (type %in% c("p", "d") & is.null(x)) {
    if (is.null(x_num)) {
      x_num <- 100
    }
    p_min <- 0
    if (object$model == "bingp") {
      p_min <- pred_pbingp(ev_obj = object, q = object$thresh,
                           n_years = n_years, npy = npy, lower_tail = TRUE)$y
    }
    ep <- c(0.001, 0.01)
    x <- rbind(pmax(ep[1], p_min), pmax(1 - ep[2], p_min))
    x <- matrix(x, nrow = 2, ncol = n_y, byrow = FALSE)
    x <- qpred(object, p = x, n_years = n_years, npy = npy,
               lower_tail = TRUE, big_q = big_q)$y
    x <- apply(x, 2, function(x) seq(from = x[1], to = x[2], len = x_num))
  }
  if (type == "p") {
    ret_obj <- ppred(object, q = x, n_years = n_years, npy = npy,
                     lower_tail = lower_tail)
  }
  if (type == "d") {
    ret_obj <- dpred(object, x = x, n_years = n_years, npy = npy,
                     log = log)
  }
  if (type == "q") {
    ret_obj <- qpred(object, p = x, n_years = n_years, npy = npy,
                     lower_tail = lower_tail, big_q = big_q)
  }
  if (type == "r") {
    ret_obj <- rpred(object, n_years = n_years, npy = npy)
  }
  ret_obj$type <- type
  ret_obj$n_years <- n_years
  ret_obj$npy <- npy
  ret_obj$level <- level
  ret_obj$hpd <- hpd
  ret_obj$lower_tail <- lower_tail
  ret_obj$log <- log
  ret_obj$data <- object$data
  ret_obj$model <- object$model
  class(ret_obj) <- "evpred"
  return(ret_obj)
}

# ----------------------------- set_npy ---------------------------------

set_npy <- function(object, npy = NULL){
  model <- object$model
  if (!is.null(object$npy) & !is.null(npy)) {
    warning(paste("Two values of npy supplied.  The value npy = ", npy,
                  " from the current call has been used.", sep=""))
  }
  if (!is.null(object$npy) & is.null(npy)) {
    npy <- object$npy
  }
  if (object$model == "bingp") {
    if (is.null(object$npy) & is.null(npy)) {
      stop("model=bingp: npy must be given, here or in call to rpost/rpost_rcpp.")
    }
  } else if (object$model %in% c("gev", "os")) {
    # If npy is not supplied and the model is GEV or OS then assume that
    # npy = 1, that is, the data were annual maxima or annual order statistics
    # respectively.
    if (is.null(npy)) {
      npy <- 1
    }
  } else if (object$model == "pp") {
    # Similarly if npy is not supplied and the model is PP then assume that
    # blocks of length one year were set by noy in the call to
    # rpost()/rpost_rcpp().
    n <- length(object$data)
    noy <- object$noy
    if (is.null(npy)) {
      npy <- n / noy
    }
  }
  return(npy)
}

# ----------------------------- dpred ---------------------------------

dpred <- function(ev_obj, x, n_years = 100, npy = NULL, log = FALSE) {
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_dgev(ev_obj = ev_obj, x = x, n_years = n_years,
                          npy = npy, log = log)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_dbingp(ev_obj = ev_obj, x = x, n_years = n_years,
                           npy = npy, log = log)
  }
  return(ret_obj)
}

# ----------------------------- ppred ---------------------------------

ppred <- function(ev_obj, q, n_years = 100, npy = NULL, lower_tail = TRUE) {
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_pgev(ev_obj = ev_obj, q = q, n_years = n_years,
                         npy = npy, lower_tail = lower_tail)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_pbingp(ev_obj = ev_obj, q = q, n_years = n_years,
                           npy = npy, lower_tail = lower_tail)
  }
  return(ret_obj)
}

# ----------------------------- qpred ---------------------------------

qpred <- function(ev_obj, p, n_years = 100, npy = NULL, lower_tail = TRUE,
                  big_q) {
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_qgev(ev_obj = ev_obj, p = p, n_years = n_years,
                         npy = npy, lower_tail = lower_tail,
                         big_q = big_q)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_qbingp(ev_obj = ev_obj, p = p, n_years = n_years,
                           npy = npy, lower_tail = lower_tail,
                           big_q = big_q)
  }
  return(ret_obj)
}

# ----------------------------- rpred ---------------------------------

rpred <- function(ev_obj, n_years = 100, npy = NULL) {
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    ret_obj <- pred_rgev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  } else if (ev_obj$model == "bingp") {
    ret_obj <- pred_rbingp(ev_obj = ev_obj, n_years = n_years, npy = npy)
  }
  return(ret_obj)
}

# ----------------------------- ipred ---------------------------------

ipred <- function(ev_obj, n_years = 100, npy = NULL, level = 95,
                  hpd = FALSE, big_q) {
  if (any(level <= 0) || any(level >= 100)) {
    stop("level must be in (0, 100)")
  }
  if (ev_obj$model %in% c("gev", "os", "pp")) {
    qfun <- pred_qgev
    pfun <- pred_pgev
    p_min <- 0
  } else if (ev_obj$model == "bingp") {
    qfun <- pred_qbingp
    pfun <- pred_pbingp
    # Find the smallest allowable value of p, i.e. the one that corresponds
    # to the threshold used in the call to rpost()/rpost_rcpp().  This enables
    # us to check whether or not that it is possible to calculate a given
    # level% predictive interval without extending below the threshold.
    p_min <- pfun(ev_obj = ev_obj, q = ev_obj$thresh, n_years = n_years,
                  npy = npy, lower_tail = TRUE)$y
  }
  # Find the equi-tailed level% interval (for hpd = FALSE)
  p1 <- (1 - level / 100) / 2.
  # Create list to return.
  ret_obj <- list()
  # If the equi-tailed interval extends below the threshold then return NAs,
  # for both long and for short.
  if (p1 < p_min) {
    ret_obj$long <- matrix(NA, ncol = 1, nrow = 2)
    ret_obj$short <- matrix(NA, ncol = 1, nrow = 2)
    return(ret_obj)
  }
  # Non-exceedance probabilities corresponding to the equi-tailed level%
  # intervals.
  pp <- c(p1, 1 - p1)
  # Find the corresponding predictive quantiles.
  ret_obj$long <- qfun(ev_obj = ev_obj, p = pp, n_years = n_years, npy = npy,
                       lower_tail = TRUE, big_q = big_q)$y
  if (!hpd) {
    ret_obj$short <- matrix(NA, ncol = 1, nrow = 2)
    return(ret_obj)
  }
  qq <- ret_obj$long
  # Set a small probability that is less than p1.
  ep <- min(.Machine$double.eps ^ 0.5, p1 / 10)
  # ... but no smaller than p_min
  ep <- max(ep, p_min)
  #
  # Assume that the hpd interval has a lower lower endpoint than the
  # equi-tailed interval.  This will typically be the case for predictive
  # intervals for N-year maxima for large N.  Then the lower endpoint of
  # the hpd interval should correspond to a non-exceedance probability
  # between ep and p1.
  #
  # Set non-exceedance probabilities for an interval that should lie below the
  # hpd interval.
  p_low <- c(ep, ep + level / 100)
  # Find the corresponding predictive quantiles.
  q_low <- qfun(ev_obj = ev_obj, p = p_low, n_years = n_years, npy = npy,
                lower_tail = TRUE, init_q = qq, big_q = big_q)$y
  # q_lower and q_upper are the respective intervals within which the lower
  # and upper limits of the hpd interval should fall.
  q_lower <- c(q_low[1], qq[1])
  q_upper <- c(q_low[2], qq[2])
  # ipred_hpd() uses monotonic cubic spline interpolation to estimate
  # the predictive intervals corresponding to given lower and upper
  # non-exceedance probabilities that are level/100 apart, and to
  # search for the shortest such interval.
  temp <- ipred_hpd(q_lower = q_lower, q_upper = q_upper, ev_obj = ev_obj,
                    n_years = n_years, level = level, npy = npy, pfun = pfun,
                    n_spline = 100)
  # Check that the returned interval is not produced by a non-exceedance
  # probability that is on the boundary of those considered.  If this is the
  # case then a matrix of NAs is returned from ipred_hpd().  If this happens
  # then we try on the other side of the boundary.
  if (is.na(temp[1, 1])) {
    p_low <- c(1 - level / 100 - ep, 1- ep)
    q_low <- qfun(ev_obj = ev_obj, p = p_low, n_years = n_years, npy = npy,
                  lower_tail = TRUE, init_q = qq)$y
    q_lower <- c(qq[1], q_low[1])
    q_upper <- c(qq[2], q_low[2])
    temp <- ipred_hpd(q_lower = q_lower, q_upper = q_upper, ev_obj = ev_obj,
                      n_years = n_years, level = level, npy = npy, pfun = pfun,
                      n_spline = 100)
  }
  ret_obj$short <- temp
  return(ret_obj)
}

ipred_hpd <- function(q_lower, q_upper, ev_obj, n_years, level, npy, pfun,
                      n_spline = 100) {
  # To find the hpd interval we want to call the quantile function
  # as little as possible because this is slow.  Therefore, we estimate the
  # quantile function in the regions that matter using monotonic cubic
  # splines fitted to pairs of fixed quantiles and the estimates of the
  # corresponding distribution function values, calculated using pfun.
  #
  # Start with sequences of quantiles that are equally-spaced.
  q_lower <- seq(q_lower[1], q_lower[2], len = n_spline)
  q_upper <- seq(q_upper[1], q_upper[2], len = n_spline)
  # Find the corresponding non-exceedance probabilities.
  p_lower <- pfun(ev_obj = ev_obj, q = q_lower, n_years = n_years,
                  npy = npy, lower_tail = TRUE)$y
  p_upper <- pfun(ev_obj = ev_obj, q = q_upper, n_years = n_years,
                  npy = npy, lower_tail = TRUE)$y
  p_range <- range(p_lower)
  # Perform monotonic cubic spline interpolation of (p,q).
  lower_spline <- stats::splinefun(x = p_lower, y = q_lower, method = "hyman")
  upper_spline <- stats::splinefun(x = p_upper, y = q_upper, method = "hyman")
  # Now set equally-spaced lower probabilities and their corresponding upper
  # probabilities for a level% interval.  This is because we will search over
  # the probabilities and the spline interpolation should work better if the
  # input probabilities are approximately equi-spaced and we can also find
  # the lengths of the level% intervals at the knots, before performing the
  # minimisation, if we do things this way.
  p_lower <- seq(p_range[1], p_range[2], len = n_spline)
  p_upper <- p_lower + level / 100
  # Use the spline interpolation to estimate the corresponding quantiles.
  lower_qs <-  lower_spline(x = p_lower)
  upper_qs <-  upper_spline(x = p_upper)
  # Find the actual (no spline) probabilities corresponding to these quantiles.
  p_lower <- pfun(ev_obj = ev_obj, q = lower_qs, n_years = n_years,
                  npy = npy, lower_tail = TRUE)$y
  p_upper <- pfun(ev_obj = ev_obj, q = upper_qs, n_years = n_years,
                  npy = npy, lower_tail = TRUE)$y
  # Find the lengths of the intervals and the location of the minimum length.
  q_length <- upper_qs - lower_qs
  where_min_p <- which.min(q_length)
  # If the best p is on the lower boundary then return these values.
  if (where_min_p == 1) {
    return(matrix(c(lower_qs[1], upper_qs[1]), ncol = 1, nrow = 2))
  }
  # If the best p is on the upper boundary then return NAs.
  if (where_min_p == n_spline) {
    return(matrix(NA, ncol = 1, nrow = 2))
  }
  # Use the value of p that gives the shortest interval as an initial estimate.
  which_p <- which.min(q_length)
  p_init <- p_lower[which_p]
  p_range <- c(p_lower[which_p - 1], p_lower[which_p + 1])
  #
  # Objective function that returns the length of the level% interval for a
  # given lower non-exceedance probability.
  ob_fun <- function(p1, ev_obj, n_years, npy, level) {
    p <- c(p1, p1 + level / 100)
    lower_limit <-  lower_spline(x = p[1])
    upper_limit <-  upper_spline(x = p[2])
    limits <- c(lower_limit, upper_limit)
    return(structure(diff(limits), limits = limits))
  }
  temp <- stats::nlminb(p_init, ob_fun, ev_obj = ev_obj, n_years = n_years,
                        npy = npy, level = level, lower = p_range[1],
                        upper = p_range[2])
  temp <- ob_fun(p1 = temp$par, ev_obj = ev_obj, n_years = n_years, npy = npy,
                 level = level)
  return(matrix(attr(temp, "limits"), ncol = 1, nrow = 2))
}

# ============================ GEV functions ============================

# ----------------------------- pred_dgev ---------------------------------

pred_dgev <- function(ev_obj, x, n_years = 100, npy = NULL, log = FALSE) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_y <- length(n_years)
  if (is.vector(x)) {
    x <- matrix(x, ncol = n_y, nrow = length(x), byrow = FALSE)
  }
  if (ncol(x) != n_y) {
    stop("quantiles must be a vector or a matrix with length(n_years) columns")
  }
  d <- x
  temp <- function(x, loc, scale, shape, m) {
    return(mean(dgev(x = x, loc = loc, scale = scale, shape = shape, m = m)))
  }
  for (i in 1:n_y) {
    # Calculate the GEV pdf at x for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
    d[, i] <- sapply(x[, i], temp, loc = loc, scale = scale, shape = shape,
                     m = mult[i])
  }
  if (log) {
    d <- log(d)
  }
  return(list(x = x, y = d))
}

# ----------------------------- pred_pgev ---------------------------------

pred_pgev <- function(ev_obj, q, n_years = 100, npy = NULL,
                      lower_tail = TRUE) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_y <- length(n_years)
  if (is.vector(q)) {
    q <- matrix(q, ncol = n_y, nrow = length(q), byrow = FALSE)
  }
  if (ncol(q) != n_y) {
    stop("quantiles must be a vector or a matrix with length(n_years) columns")
  }
  p <- q
  temp <- function(q, loc, scale, shape, m) {
    return(mean(pgev(q = q, loc = loc, scale = scale, shape = shape, m = m)))
  }
  for (i in 1:n_y) {
    # Calculate the GEV cdf at q for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
    p[, i] <- sapply(q[, i], temp, loc = loc, scale = scale, shape = shape,
                     m = mult[i])
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  return(list(x = q, y = p))
}

# ----------------------------- pred_qgev ---------------------------------

pred_qgev <- function(ev_obj, p, n_years = 100, npy = NULL,
                      lower_tail = TRUE, init_q = NULL, big_q) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  if (!lower_tail) {
    p <- 1 - p
  }
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_y <- length(n_years)
  if (is.vector(p)) {
    p <- matrix(p, ncol = n_y, nrow = length(p), byrow = FALSE)
  }
  if (ncol(p) != n_y) {
    stop("quantiles must be a vector or a matrix with length(n_years) columns")
  }
  n_p <- nrow(p)
  q <- p
  # Check that the dimensions of init_q are OK.
  ok_init_q <- FALSE
  if (!is.null(init_q)) {
    ok_init_q <- init_q_check(init_q = init_q, n_p = n_p, n_y = n_y)
  }
  if (!ok_init_q) {
    init_q <- matrix(NA, nrow = n_p, ncol = n_y)
    temp <- function(p, loc, scale, shape, m) {
      return(mean(qgev(p = p, loc = loc, scale = scale, shape = shape, m = m)))
    }
  }
  lower <- min(ev_obj$data)
  upper <- big_q
  u_minus_l <- upper - lower
  for (i in 1:n_y) {
    # Calculate the GEV quantile at p for each combination of (loc, scale, shape)
    # in the posterior sample, and take the mean.
    #
    # This gives reasonable initial estimates for the predictive quantiles.
    if (!ok_init_q) {
      init_q[, i] <- sapply(p[, i], temp, loc = loc, scale = scale,
                            shape = shape, m = mult[i])
    }
    #
    ob_fn <- function(q, ev_obj, p, n_years, npy) {
      p_val <- pred_pgev(ev_obj = ev_obj, q = q, n_years = n_years,
                         npy = npy)$y
      return(p_val - p)
    }
    for (j in 1:n_p) {
      f_upper <- ob_fn(upper, ev_obj = ev_obj, p = p[j, i],
                       n_years = n_years[i], npy = npy)
      k <- 1
      while (f_upper < 0) {
        upper <- lower + u_minus_l * (10 ^ k)
        k <- k + 1
        f_upper <- ob_fn(upper, ev_obj = ev_obj, p = p[j, i],
                         n_years = n_years[i], npy = npy)
      }
      qtemp <- stats::uniroot(f = ob_fn, ev_obj = ev_obj, p = p[j, i],
                              n_years = n_years[i], npy = npy,
                              lower = lower, upper = upper, f.upper = f_upper,
                              tol = .Machine$double.eps^0.5)
      q[j, i] <- qtemp$root
    }
  }
  return(list(x = p, y = q))
}

# ----------------------------- init_q_check ---------------------------------

init_q_check <- function(init_q, n_p, n_y) {
  ok_init_q <- TRUE
  if (is.vector(init_q)) {
    if (length(init_q) == 1) {
      q_mat <- matrix(init_q, ncol = 1, nrow = n_p)
    } else if (length(init_q) == n_p) {
      q_mat <- matrix(init_q)
    } else {
      warning("init_q has an invalid size: initial values set internally.")
      ok_init_q <- FALSE
    }
  } else if (is.matrix(init_q)) {
    if (nrow(init_q) == n_p) {
      q_mat <- matrix(init_q, ncol = n_y, nrow = n_p, byrow = FALSE)
    } else {
      warning("dim(init_q) is invalid: initial values set internally.")
      ok_init_q <- FALSE
    }
  } else {
    warning("init_q has an invalid type: initial values set internally.")
    ok_init_q <- FALSE
  }
  return(ok_init_q)
}

# ----------------------------- pred_rgev ---------------------------------

pred_rgev <- function(ev_obj, n_years = 100, npy = NULL) {
  # Determine the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  mult <- setup_pred_gev(ev_obj = ev_obj, n_years = n_years, npy = npy)
  #
  loc <- ev_obj$sim_vals[, 1]
  scale <- ev_obj$sim_vals[, 2]
  shape <- ev_obj$sim_vals[, 3]
  n_sim <- length(loc)
  n_y <- length(n_years)
  r_mat <- matrix(NA, nrow = n_sim, ncol = n_y)
  for (i in 1:n_y) {
    # Simulate a single observation from a GEV distribution corresponding
    # to each parameter combination in the posterior sample.
    r_mat[, i] <- rgev(n = n_sim, loc = loc, scale = scale, shape = shape,
                       m = mult[i])
  }
  return(list(y = r_mat))
}

# -------------------------- setup_pred_gev ------------------------------

setup_pred_gev <- function(ev_obj, n_years, npy) {
  #
  # Determines the number, mult, of blocks in n_years years, so that
  # the GEV parameters can be converted to the n_years level of aggregation.
  #
  # Args:
  #   ev_obj  : Object of class evpost return by rpost()/rpost_rcpp().
  #   n_years : A numeric vector. Values of N.
  #   npy     : The mean number of observations per year of data, after
  #             excluding any missing values.  npy has either been supplied
  #             by the user or, if this is not the case, was set using
  #             set_npy() based on assumptions that in the original call to
  #             rpost()/rpost_rcpp() the GEV parameters relate to blocks
  #             of length one year.
  #
  # Returns: A numeric scalar.  The value of mult.
  #
  if (ev_obj$model %in% c("gev", "os")) {
    mult <- n_years * npy
  }
  if (ev_obj$model == "pp") {
    n <- length(ev_obj$data)
    noy <- ev_obj$noy
    mult <- n_years * npy * noy / n
  }
  return(mult)
}

# ============================ binGP functions ============================

# ----------------------------- pred_dbingp ---------------------------------

pred_dbingp <- function(ev_obj, x, n_years = 100, npy = NULL,
                        log = FALSE) {
  # Check that q is not less than the threshold used in the call to
  # rpost()/rpost_rcpp().
  thresh <- ev_obj$thresh
  if (any(x < thresh)) {
    stop("Invalid x: no element of x can be less than the threshold.")
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  n_y <- length(n_years)
  if (is.vector(x)) {
    x <- matrix(x, ncol = n_y, nrow = length(x), byrow = FALSE)
  }
  if (ncol(x) != n_y) {
    stop("quantiles must be a vector or a matrix with length(n_years) columns")
  }
  d <- x
  temp <- function(x, p_u, scale, shape, thresh, mult) {
    # Calculate the distribution function of raw observations, evaluated at q.
    raw_df <- pbingp(q = x, p_u = p_u, loc = thresh, scale = scale,
                     shape = shape)
    # Evaluate the derivative of raw_df ^ mult with respect to x.
    t1 <- mult * exp((mult - 1) * log(raw_df))
    t2 <- dbingp(x = x, p_u = p_u, loc = thresh, scale = scale, shape = shape)
    # Return the mean of the posterior sample.
    return(mean(t1 * t2))
  }
  # For each value in n_years calculate the distribution function of the
  # n_years maximum.
  mult <- npy * n_years
  for (i in 1:n_y) {
    d[, i] <- sapply(x[, i], temp, p_u = p_u, scale = scale, shape = shape,
                     thresh = thresh, mult = mult[i])
  }
  if (log) {
    d <- log(d)
  }
  return(list(x = x, y = d))
}

# ----------------------------- pred_pbingp ---------------------------------

pred_pbingp <- function(ev_obj, q, n_years = 100, npy = NULL,
                        lower_tail = TRUE) {
  # Check that q is not less than the threshold used in the call to
  # rpost()/rpost_rcpp().
  thresh <- ev_obj$thresh
  if (any(q < thresh)) {
    stop("Invalid q: no element of q can be less than the threshold.")
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  n_y <- length(n_years)
  if (is.vector(q)) {
    q <- matrix(q, ncol = n_y, nrow = length(q), byrow = FALSE)
  }
  if (ncol(q) != n_y) {
    stop("quantiles must be a vector or a matrix with length(n_years) columns")
  }
  p <- q
  temp <- function(q, p_u, scale, shape, thresh, mult) {
    # Calculate the distribution function of raw observations, evaluated at q.
    raw_df <- pbingp(q = q, p_u = p_u, loc = thresh, scale = scale,
                      shape = shape)
    # Raise this to the power of mult to find the distribution function of
    # the n_year maximum.
    return(mean(exp(mult * log(raw_df))))
  }
  # For each value in n_years calculate the distribution function of the
  # n_years maximum.
  mult <- npy * n_years
  for (i in 1:n_y) {
    p[, i] <- sapply(q[, i], temp, p_u = p_u, scale = scale, shape = shape,
                     thresh = thresh, mult = mult[i])
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  return(list(x = q, y = p))
}

# ----------------------------- pred_qbingp ---------------------------------

pred_qbingp <- function(ev_obj, p, n_years = 100, npy = NULL,
                        lower_tail = TRUE, init_q = NULL, big_q) {
  if (!lower_tail) {
    p <- 1 - p
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  n_y <- length(n_years)
  if (is.vector(p)) {
    p <- matrix(p, ncol = n_y, nrow = length(p), byrow = FALSE)
  }
  if (ncol(p) != n_y) {
    stop("quantiles must be a vector or a matrix with length(n_years) columns")
  }
  # Check that p is not less than the binGP predictive distribution function
  # evaluated at the threshold used in the call to rpost()/rpost_rcpp().
  thresh <- ev_obj$thresh
  p_ok <- rep(TRUE, n_y)
  for (i in 1:n_y) {
    p_min <- pred_pbingp(ev_obj = ev_obj, q = thresh, n_years = n_years[i],
                         npy = npy, lower_tail = lower_tail)$y
    p_min <- as.vector(p_min)
    if (any(p[, i] < p_min)) {
      p_ok[i] <- FALSE
    }
  }
  if (any(!p_ok)) {
    stop("Invalid p for at least one value in n_years")
  }
  n_p <- nrow(p)
  q <- p
  # Check that the dimensions of init_q are OK.
  ok_init_q <- FALSE
  if (!is.null(init_q)) {
    ok_init_q <- init_q_check(init_q = init_q, n_p = n_p, n_y = n_y)
  }
  if (!ok_init_q) {
    init_q <- matrix(NA, nrow = n_p, ncol = n_y)
    temp <- function(p, p_u, loc, scale, shape, mult) {
      pnew <- exp(log(p) / mult)
      # We need to avoid calling qbingp with a probability that is lower than
      # 1 - p_u, i.e. corresponds to a quantile that this below the threshold.
      co <- pnew < 1 - p_u
      qv <- qbingp(pnew, p_u = p_u[!co], loc = loc,
                        scale = scale[!co], shape = shape[!co])
      return(mean(qv))
    }
  }
  mult <- npy * n_years
  lower <- ev_obj$thresh
  upper <- big_q
  u_minus_l <- upper - lower
  for (i in 1:n_y) {
    # Calculate the quantile at p for each combination of parameters in the
    # in the posterior sample, and take the mean.
    #
    # This gives reasonable initial estimates for the predictive quantiles.
    if (!ok_init_q) {
      init_q[, i] <- sapply(p[, i], temp, p_u = p_u, loc = thresh,
                            scale = scale, shape = shape, mult = mult[i])
    }
    #
    ob_fn <- function(q, ev_obj, p, n_years, npy) {
      p_val <- pred_pbingp(ev_obj = ev_obj, q = q, n_years = n_years,
                           npy = npy)$y
      return(p_val - p)
    }
    for (j in 1:n_p) {
      f_upper <- ob_fn(upper, ev_obj = ev_obj, p = p[j, i],
                       n_years = n_years[i], npy = npy)
      k <- 1
      while (f_upper < 0) {
        upper <- lower + u_minus_l * (10 ^ k)
        k <- k + 1
        f_upper <- ob_fn(upper, ev_obj = ev_obj, p = p[j, i],
                         n_years = n_years[i], npy = npy)
      }
      # Note: pred_pbingp() cannot be evaluated for q < ev_obj$thresh.
      qtemp <- stats::uniroot(f = ob_fn, ev_obj = ev_obj, p = p[j, i],
                              n_years = n_years[i], npy = npy,
                              lower = lower, upper = upper, f.upper = f_upper,
                              tol = .Machine$double.eps^0.5)
      q[j, i] <- qtemp$root
    }
  }
  return(list(x = p, y = q))
}

# ----------------------------- pred_rbingp ---------------------------------

pred_rbingp <- function(ev_obj = ev_obj, n_years = n_years, npy = npy) {
  if (!inherits(ev_obj, "evpost")) {
    stop("ev_obj must be an object produced by rpost() or rpost_rcpp()")
  }
  if (ev_obj$model != "bingp") {
    stop("The model in the call to rpost() or rpost_rcpp() must be bingp.")
  }
  if (is.null(npy)) {
    stop("npy must be supplied.")
  }
  # For each value in n_years find the value p_min below which the
  # correspnding simulated value is below the threshold.  A missing value
  # NA is returned to indicate this.
  n_y <- length(n_years)
  thresh <- ev_obj$thresh
  for (i in 1:n_y) {
    p_min <- pred_pbingp(ev_obj = ev_obj, q = thresh, n_years = n_years[i],
                         npy = npy, lower_tail = TRUE)$y
    p_min <- as.vector(p_min)
  }
  # Extract posterior sample of parameters p_u, sigma_u, xi.
  p_u <- ev_obj$bin_sim_vals
  scale <- ev_obj$sim_vals[, 1]
  shape <- ev_obj$sim_vals[, 2]
  # Set up a matrix to contain the results.
  n_sim <- length(p_u)
  r_mat <- matrix(NA, nrow = n_sim, ncol = n_y)
  #
  mult <- npy * n_years
  # We use the sample underlying simulated uniforms for each value in n_years.
  u <- stats::runif(n_sim)
  for (i in 1:n_y) {
    x_val <- (1 - u ^ (1 / mult[i])) / p_u
    new_mult <- box_cox_vec(x = x_val, lambda = -shape)
    r_mat[, i] <- ifelse(u < p_min, NA, thresh - scale * new_mult)
  }
  return(list(y = r_mat))
}
