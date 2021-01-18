# ================================ set_prior ===============================

#' Construction of prior distributions for extreme value model parameters
#'
#' Constructs a prior distribution for use as the argument \code{prior} in
#' \code{\link{rpost}} and \code{\link{rpost_rcpp}}.  The user can either
#' specify their own prior function, returning the \strong{log} of the prior
#' density, (using an R function or an external pointer to a compiled C++
#' function) and arguments for hyperparameters or choose from a list of
#' in-built model-specific priors.  Note that the arguments
#' \code{model = "gev"}, \code{model = "pp"} and \code{model =="os"} are
#' equivalent because a prior is specified is the GEV parameterisation in each
#' of these cases.
#' Note also that for \code{model = "pp"} the prior GEV parameterisation
#' relates to the value of \code{noy} subsequently supplied to
#' \code{\link{rpost}} or \code{\link{rpost_rcpp}}.
#' The argument \code{model} is used for consistency with \code{rpost}.
#'
#' @param prior Either
#' \itemize{
#'   \item {An R function, or a pointer to a user-supplied compiled
#'   C++ function, that returns the value of the log of the prior density
#'   (see \strong{Examples}), or}
#'   \item {A character string giving the name of the prior.
#'     See \strong{Details} for a list of priors available for each model.}
#' }
#' @param model A character string.  If \code{prior} is a character string
#'   then \code{model} gives the extreme value model to be used.  Using
#'   either \code{model = "gev"}, \code{model = "pp"} or
#'   \code{model = "os"} will result in the same (GEV) parameterisation.
#'   If \code{prior} is a function then the value of \code{model} is stored
#'   so that in the subsequent call to \code{rpost}, consistency of the
#'   prior and extreme value model parameterisations can be checked.
#' @param ... Further arguments to be passed to the user-supplied or
#'   in-built prior function.  For details of the latter see \strong{Details}
#'   and/or the relevant underlying function: \code{\link{gp_norm}},
#'   \code{\link{gp_mdi}}, \code{\link{gp_flat}}, \code{\link{gp_flatflat}},
#'   \code{\link{gp_jeffreys}}, \code{\link{gp_beta}},
#'   \code{\link{gev_norm}}, \code{\link{gev_loglognorm}},
#'   \code{\link{gev_mdi}}, \code{\link{gev_flat}}, \code{\link{gev_flatflat}},
#'   \code{\link{gev_beta}}, \code{\link{gev_prob}}, \code{\link{gev_quant}}.
#'   All these priors have the arguments \code{min_xi} (prior lower bound on
#'   \eqn{\xi}) and \code{max_xi} (prior upper bound on \eqn{\xi}).
#' @details Of the in-built named priors available in revdbayes only
#'   those specified using \code{prior = "prob"} (\code{\link{gev_prob}}),
#'   \code{prior = "quant"} (\code{\link{gev_quant}})
#'   \code{prior = "norm"} (\code{\link{gev_norm}}) or
#'   \code{prior = "loglognorm"} (\code{\link{gev_loglognorm}}) are proper.
#'   If \code{model = "gev"} these priors are equivalent to priors available
#'   in the evdbayes package, namely \code{\link[evdbayes:prior]{prior.prob}},
#'   \code{\link[evdbayes:prior]{prior.quant}},
#'   \code{\link[evdbayes:prior]{prior.norm}} and
#'   \code{\link[evdbayes:prior]{prior.loglognorm}}.
#'
#'   The other in-built priors are improper, that is, the integral of the
#'   prior function over its support is not finite.  Such priors do not
#'   necessarily result in a proper posterior distribution. Northrop and
#'   Attalides (2016) consider the issue of posterior propriety in Bayesian
#'   extreme value analyses.  In most of improper priors below the prior for
#'   the scale parameter \eqn{\sigma} is taken to be \eqn{1/\sigma},
#'   i.e. a flat prior for \eqn{\log \sigma}{log \sigma}.  Here we denote the scale
#'   parameter of the GP distribution by \eqn{\sigma}, whereas we use
#'   \eqn{\sigma_u} in the revdbayes vignette.
#'
#'   For all in-built priors the arguments \code{min_xi} and \code{max_xi} may
#'   be supplied by the user.  The prior density is set to zero for any value
#'   of the shape parameter \eqn{\xi} that is outside
#'   (\code{min_xi}, \code{max_xi}).  This will override the default values
#'   of \code{min_xi} and \code{max_xi} in the named priors detailed above.
#'
#'   \strong{Extreme value priors.} It is typical to use either
#'   \code{prior = "prob"} (\code{\link{gev_prob}}) or
#'   \code{prior = "quant"} (\code{\link{gev_quant}}) to set an informative
#'   prior and one of the other prior (or a user-supplied function) otherwise.
#'   The names of the in-built extreme value priors set using \code{prior}
#'   and details of hyperparameters are:
#' \itemize{
#'   \item{\code{"prob"}.  A prior for GEV parameters \eqn{(\mu, \sigma, \xi)}
#'     based on Crowder (1992).  See \code{\link{gev_prob}} for details.
#'     See also Northrop et al. (2017) and Stephenson (2016).
#'   }
#'   \item{\code{"quant"}.  A prior for GEV parameters \eqn{(\mu, \sigma, \xi)}
#'     based on Coles and Tawn (1996). See \code{\link{gev_quant}} for details.
#'   }
#'   \item {\code{"norm"}.
#'
#'   For \code{model = "gp"}:
#'     (\eqn{\log \sigma, \xi}{log \sigma, \xi}), is bivariate normal
#'     with mean \code{mean} (a numeric vector of length 2) and covariance
#'     matrix \code{cov} (a symmetric positive definite 2 by 2 matrix).
#'
#'   For \code{model = "gev"}:
#'     (\eqn{\mu, \log \sigma, \xi}{\mu, log \sigma, \xi}), is trivariate
#'     normal with mean \code{mean} (a numeric vector of length 3) and
#'     covariance matrix \code{cov} (a symmetric positive definite 3 by 3
#'     matrix).
#'   }
#'   \item {\code{"loglognorm"}.  For \code{model = "gev"} only:
#'     (\eqn{\log \mu, \log \sigma, \xi}{log \mu, log \sigma, \xi}), is
#'     trivariate normal with mean \code{mean} (a numeric vector of length 3)
#'     and covariance matrix \code{cov} (a symmetric positive definite 3 by 3
#'     matrix).
#'   }
#'   \item {\code{"mdi"}.
#'
#'   For \code{model = "gp"}: (an extended version
#'     of) the maximal data information (MDI) prior, that is,
#'     \deqn{\pi(\sigma, \xi) = \sigma^{-1} \exp[-a(\xi + 1)], {\rm ~for~}
#'     \sigma > 0, \xi \geq -1, a \geq 0.}{%
#'     \pi(\sigma, \xi) = (1/ \sigma) exp[- a (\xi + 1)], for
#'     \sigma > 0, \xi >= -1, a >= 0.}
#'     The value of \eqn{a} is set using the argument \code{a}.  The default
#'     value is \eqn{a = 1}, which gives the MDI prior.
#'
#'     For \code{model = "gev"}: (an extended version
#'     of) the maximal data information (MDI) prior, that is,
#'     \deqn{\pi(\mu, \sigma, \xi) = \sigma^{-1} \exp[-a(\xi + 1)],
#'     {\rm ~for~} \sigma > 0, \xi \geq -1, a \geq 0.}{%
#'     \pi(\mu, \sigma, \xi) = (1/ \sigma) exp[- a (\xi + 1)], for
#'     \sigma > 0, \xi >= -1, a >= 0.}
#'     The value of \eqn{a} is set using the argument \code{a}.  The default
#'     value is \eqn{a = \gamma}, where \eqn{\gamma = 0.57721} is Euler's
#'     constant, which gives the MDI prior.
#'
#'     For each of these cases \eqn{\xi} must be is bounded below
#'     \emph{a priori} for the posterior to be proper
#'     (Northrop and Attalides, 2016).  An argument for the
#'     bound \eqn{\xi \geq -1}{\xi >= -1} is that for \eqn{\xi < -1} the
#'     GP (GEV) likelihood is unbounded above as \eqn{-\sigma/\xi}
#'     (\eqn{\mu - \sigma/\xi})) approaches the sample maximum.  In
#'     maximum likelihood estimation of GP parameters (Grimshaw, 1993)
#'     and GEV parameters a local maximum of the likelihood
#'     is sought on the region
#'     \eqn{\sigma > 0, \xi \geq -1}{\sigma > 0, \xi >= -1}.
#'   }
#'   \item{\code{"flat"}.
#'
#'     For \code{model = "gp"}: a flat prior for
#'     \eqn{\xi} (and for \eqn{\log \sigma}{log \sigma}):
#'     \deqn{\pi(\sigma, \xi) = \sigma^{-1}, {\rm ~for~} \sigma > 0.}{%
#'           \pi(\sigma, \xi) = (1/ \sigma), for \sigma > 0.}
#'
#'     For \code{model = "gev"}: a flat prior for
#'     \eqn{\xi} (and for \eqn{\mu} and \eqn{\log \sigma}{log \sigma}):
#'     \deqn{\pi(\mu, \sigma, \xi) = \sigma^{-1}, {\rm ~for~} \sigma > 0.}{%
#'           \pi(\mu, \sigma, \xi) = (1/ \sigma), for \sigma > 0.}
#'   }
#'   \item{\code{"flatflat"}.
#'
#'     For \code{model = "gp"}: flat priors for
#'     \eqn{\sigma} and \eqn{\xi}:
#'     \deqn{\pi(\sigma, \xi) = 1, {\rm ~for~} \sigma > 0.}{%
#'           \pi(\sigma, \xi) = 1, for \sigma > 0.}
#'
#'     For \code{model = "gev"}: flat priors for \eqn{\mu}, \eqn{\sigma}
#'     and \eqn{\xi}:
#'     \deqn{\pi(\mu, \sigma, \xi) = 1, {\rm ~for~} \sigma > 0.}{%
#'           \pi(\mu, \sigma, \xi) = 1, for \sigma > 0.}
#'
#'     Therefore, the posterior is proportional to the likelihood.
#'   }
#'   \item{\code{"jeffreys"}.  For \code{model = "gp"} only: the Jeffreys
#'     prior (Castellanos and Cabras, 2007):
#'     \deqn{\pi(\sigma, \xi) = \sigma^{-1}(1+\xi)^{-1}(1+2\xi)^{-1/2},
#'       {\rm ~for~} \sigma > 0, \xi > -1 / 2.}{%
#'       \pi(\sigma, \xi) = 1/ [\sigma (1+\xi) \sqrt(1+2\xi)],
#'       for \sigma > 0, \xi > -1 / 2.}
#'
#'     In the GEV case the Jeffreys prior doesn't yield a proper posterior
#'     for any sample size.  See Northrop and Attalides (2016) for details.
#'   }
#'   \item{\code{"beta"}.
#'     For \code{model = "gp"}: a beta-type(p, q)
#'     prior is used for xi on the interval (\code{min_xi}, \code{max_xi}):
#'     \deqn{\pi(\sigma, \xi) = \sigma^{-1} (\xi - {\min}_{\xi}) ^ {p-1}
#'           ({\max}_{\xi} - \xi) ^ {q-1}, {\rm ~for~}
#'           {\min}_{\xi} < \xi < {\max}_{\xi}.}{%
#'           \pi(\sigma, \xi) = (1/\sigma) (\xi - min_xi) ^ (p-1)
#'           (max_xi - \xi) ^ (q-1), for min_xi < \xi < max_xi.}
#'
#'     For \code{model = "gev"}: similarly ...
#'     \deqn{\pi(\mu, \sigma, \xi) = \sigma^{-1} (\xi - {\min}_{\xi}) ^ {p-1}
#'           ({\max}_{\xi} - \xi) ^ {q-1}, {\rm ~for~}
#'           {\min}_{\xi} < \xi < {\max}_{\xi}.}{%
#'           \pi(\mu, \sigma, \xi) = (1/\sigma) (\xi - min_xi) ^ (p-1)
#'           (max_xi - \xi) ^ (q-1), for min_xi < \xi < max_xi.}
#'
#'     The argument \code{pq} is a vector containing \code{c(p,q)}.
#'     The default settings for this prior are \code{p = 6, q = 9} and
#'     \code{min_xi = -1/2, max_xi = 1/2}, which corresponds to the
#'     prior for \eqn{\xi} proposed in Martins and Stedinger (2000, 2001).
#'   }
#' }
#' @return A list with class \code{"evprior"}.  The first component is the
#'   input prior, i.e. either the name of the prior or a user-supplied
#'   function.  The remaining components contain the numeric values of any
#'   hyperparameters in the prior.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @seealso \code{\link{create_prior_xptr}} for creating an external
#'   pointer to a C++ function to evaluate the log-prior density.
#' @seealso \code{\link{rprior_prob}} and \code{\link{rprior_quant}} for
#'   sampling from informative prior distributions for GEV parameters.
#' @seealso \code{\link{gp_norm}}, \code{\link{gp_mdi}},
#'   \code{\link{gp_flat}}, \code{\link{gp_flatflat}},
#'   \code{\link{gp_jeffreys}}, \code{\link{gp_beta}} to see the arguments
#'   for priors for GP parameters.
#' @seealso \code{\link{gev_norm}}, \code{\link{gev_loglognorm}},
#'   \code{\link{gev_mdi}}, \code{\link{gev_flat}}, \code{\link{gev_flatflat}},
#'   \code{\link{gev_beta}}, \code{\link{gev_prob}}, \code{\link{gev_quant}}
#'   to see the arguments for priors for GEV parameters.
#' @seealso \code{\link[evdbayes:prior]{prior.prob}},
#'   \code{\link[evdbayes:prior]{prior.quant}},
#'   \code{\link[evdbayes:prior]{prior.norm}}
#'   and \code{\link[evdbayes:prior]{prior.loglognorm}} for setting a prior
#'   distribution using the evdbayes package.
#' @seealso \code{\link[evdbayes]{posterior}} for sampling from an extreme
#'   value posterior using the evdbayes package.
#' @references Castellanos, E. M. and Cabras, S. (2007) A default Bayesian
#'   procedure for the generalized Pareto distribution.
#'   \emph{Journal of Statistical Planning and Inference} \strong{137(2)},
#'   473-483. \url{https://doi.org/10.1016/j.jspi.2006.01.006}.
#' @references Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of
#'   extreme rainfall data. \emph{Appl. Statist.}, \strong{45}, 463-478.
#' @references Crowder, M. (1992) Bayesian priors based on parameter
#'   transformation using the distribution function
#'   \emph{Ann. Inst. Statist. Math.}, \strong{44}, 405-416.
#'   \url{https://link.springer.com/article/10.1007/BF00050695}.
#' @references Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
#'   for the Generalized Pareto Distribution.  \emph{Technometrics},
#'   \strong{35(2)}, 185-191.
#'   \url{https://doi.org/10.1080/00401706.1993.10485040}.
#' @references Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and
#'   Quantile Estimation for the Generalized Pareto Distribution.
#'   \emph{Technometrics}, \strong{29(3)}, 339-349.
#'   \url{https://doi.org/10.2307/1269343}.
#' @references Martins, E. S. and Stedinger, J. R. (2000) Generalized maximum
#'   likelihood generalized extreme value quantile estimators for hydrologic
#'   data, \emph{Water Resources Research}, \strong{36(3)}, 737-744.
#'   \url{https://doi.org/10.1029/1999WR900330}.
#' @references Martins, E. S. and Stedinger, J. R. (2001) Generalized maximum
#'   likelihood Pareto-Poisson estimators for partial duration series,
#'   \emph{Water Resources Research}, \strong{37(10)}, 2551-2557.
#'   \url{https://doi.org/10.1029/2001WR000367}.
#' @references Northrop, P.J. and Attalides, N. (2016) Posterior propriety in
#'   Bayesian extreme value analyses using reference priors
#'   \emph{Statistica Sinica}, \strong{26}(2), 721--743
#'   \url{https://doi.org/10.5705/ss.2014.034}.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
#' @references Stephenson, A. (2016) Bayesian inference for extreme value
#'   modelling.  In \emph{Extreme Value Modeling and Risk Analysis: Methods
#'   and Applications} (eds D. K. Dey and J. Yan), 257-280, Chapman and Hall,
#'   London. \url{https://doi.org/10.1201/b19721}.
#' @examples
#' # Normal prior for GEV parameters (mu, log(sigma), xi).
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' pn
#'
#' # Prior for GP parameters with flat prior for xi on (-1, infinity).
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' fp
#'
#' # A user-defined prior (see the vignette for details).
#' u_prior_fn <- function(x, ab) {
#'   if (x[1] <= 0 | x[2] <= -1 | x[2] >= 1) {
#'     return(-Inf)
#'   }
#'   return(-log(x[1]) + (ab[1] - 1) * log(1 + x[2]) +
#'          (ab[2] - 1) * log(1 - x[2]))
#' }
#' up <- set_prior(prior = u_prior_fn, ab = c(2, 2), model = "gp")
#'
#' # A user-defined prior using a pointer to a C++ function
#' ptr_gp_flat <- create_prior_xptr("gp_flat")
#' u_prior_ptr <- set_prior(prior = ptr_gp_flat, model = "gp")
#' @export
set_prior <- function(prior = c("norm", "loglognorm", "mdi", "flat",
                                "flatflat", "jeffreys", "beta", "prob",
                                "quant"),
                      model = c("gev", "gp", "pp", "os"), ...) {
  if (length(model) > 1) {
    warning("model not supplied: model == \"gev\" has been used.",
            immediate. = TRUE)
  }
  model <- match.arg(model)
  # If prior is a function then just return it in the required format,
  # together with any additional arguments from ....
  if (is.function(prior)) {
    temp <- list(prior = prior, ...)
    # Add trendsd to the prior list so that the prior will work with the
    # evdbayes function posterior().
    temp$trendsd <- 0
    # For the same reason we create a new prior function that has
    # trendsd as an argument.
    new_prior <- function(pars, ..., trendsd = 0) {
      return(prior(pars, ...))
    }
    temp$prior <- new_prior
    return(structure(temp, class = "evprior", model = model))
  }
  # If prior is a pointer to an external C++ function then just return it in
  # the required format.  We don't need to worry about adding the argument
  # trendsd because we don't need the prior to work with the evdbayes package.
  if (class(prior) == "externalptr") {
    temp <- list(prior = prior, ...)
    return(structure(temp, class = "evprior", model = model))
  }
  # Otherwise, call the appropriate function to set the prior with name prior.
  prior <- match.arg(prior)
  temp <- switch(model, gp = gp_prior(prior, ...), gev = gev_prior(prior, ...),
          pp = gev_prior(prior, ...), os = gev_prior(prior, ...))
  return(structure(temp, class = "evprior", model = model))
}

# ================================= GP priors =================================

gp_prior <- function(prior = c("norm", "mdi", "flat", "flatflat", "jeffreys",
                               "beta"), ...) {
  prior <- match.arg(prior)
  temp <- list(prior = paste("gp_", prior, sep=""), ...)
  # For v1.2.0 ...
  # Add default value for min_xi in Jeffreys prior if it isn't present.
  if (prior == "jeffreys" & is.null(temp$min_xi)) {
    temp$min_xi = -1/2
  }
  # Add default value for min_xi in the MDI prior if it isn't present.
  if (prior == "mdi" & is.null(temp$min_xi)) {
    temp$min_xi = -1
  }
  # ... and similarly for beta prior.
  if (prior == "beta" & is.null(temp$min_xi)) {
    temp$min_xi = -1/2
  }
  if (prior == "beta" & is.null(temp$max_xi)) {
    temp$max_xi = 1/2
  }
  # Force min_xi and max_xi to be present in the list temp.
  temp$min_xi <- max(temp$min_xi, -Inf)
  temp$max_xi <- min(temp$max_xi, Inf)
  # ... add default value for a in MDI prior if it isn't present.
  if (prior == "mdi" & is.null(temp$a)) {
    temp$a <- 1
  }
  # ... add default value for pq in beta-type prior if it isn't present.
  if (prior == "beta" & is.null(temp$pq)) {
    temp$pq = c(6, 9)
  }
  # Check for unused hyperparameter names and drop them
  hpar_vec <- switch(prior, norm = c("mean", "cov"), mdi = "a",
                     flat = NULL, jeffreys = NULL, beta = "pq")
  hpar_vec <- c(hpar_vec, "min_xi", "max_xi", "upper")
  temp <- hpar_drop(temp, hpar_vec)
  # Check for problems with min_xi and/or max_xi
  if (!is.null(temp$min_xi) & !is.null(temp$max_xi)) {
    if (temp$min_xi >= temp$max_xi)
        stop("min_xi must be less than max_xi")
  }
  if (!is.null(temp$min_xi)) {
    if (prior == "mdi" & is.infinite(temp$min_xi)) {
      stop("If min_xi=-Inf then the MDI posterior is not proper: min_xi must
             be finite.")
    }
    if (prior == "jeffreys" & temp$min_xi < -1 / 2 ) {
      temp$min_xi <- -1 / 2
      warning("min_xi < -1/2 does not make sense for the Jeffreys' prior.
               min_xi = -1/2 has been used.")
    }
  }
  # Check admissibility of hyperparameters
  if (prior == "norm") {
    mean <- temp$mean
    cov <- temp$cov
    if (length(mean) != 2 | mode(mean) != "numeric")
        stop("mean must be a numeric vector of length two")
    if (is.null(cov))
      stop("cov must be supplied")
    if (!is.matrix(cov) | any(dim(cov) != 2) | mode(cov) != "numeric")
        stop("cov must be a symmetric two by two matrix")
    if (any(abs(cov - t(cov)) > .Machine$double.eps ^ 0.5))
        warning("cov may not be symmetric")
    eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
    if (any(eg <= 0))
        warning("cov may not be positive definite")
    icov <- solve(cov)
    icov <- icov[row(icov) >= col(icov)]
    temp$cov <- NULL
    temp <- c(temp, list(icov=icov))
  }
  if (prior == "mdi" & !is.null(temp$a)) {
    a <- temp$a
    if (length(a) != 1 | !is.numeric(a) | a < 0)
        stop("a must be a non-negative numeric vector of length 1")
  }
  if (prior == "beta" & !is.null(temp$pq)) {
    pq <- temp$pq
    if (length(pq) != 2 | !is.numeric(pq) | any(pq <= 0) )
        stop("pq must be a non-negative numeric vector of length 2")
  }
  # Add trendsd to the prior list so that the prior will work with the
  # evdbayes function posterior().
  temp$trendsd <- 0
  return(temp)
}

# ------------------------------ specific GP priors -------------------------- #

#' Bivariate normal prior for GP parameters (\eqn{log \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 2.
#'   GP parameters (\eqn{\sigma, \xi}).
#' @param mean A numeric vector of length 2.  Prior mean.
#' @param icov A 2x2 numeric matrix.
#'   The inverse of the prior covariance matrix.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gp_norm <- function(pars, mean, icov, min_xi = -Inf, max_xi = Inf,
                    trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  pars[1] <- log(pars[1])
  cpar <- pars - mean
  ld <- icov[1] * cpar[1] ^ 2 + 2 * icov[2] * cpar[1] * cpar[2] +
    icov[3] * cpar[2] ^ 2
  return(-ld / 2 - pars[1])
}

#' Maximal data information (MDI) prior for GP parameters
#' (\eqn{\sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GP parameters (\eqn{\sigma, \xi}).
#' @param a A numeric scalar.  The default value, Euler's constant, gives the
#'   MDI prior.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#'   See Northrop and Attalides (2016) for details.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @references Northrop, P.J. and Attalides, N. (2016) Posterior propriety in
#'   Bayesian extreme value analyses using reference priors
#'   \emph{Statistica Sinica}, \strong{26(2)}, 721--743
#'   \url{https://doi.org/10.5705/ss.2014.034}.
#' @export
gp_mdi <- function(pars, a = 1, min_xi = -1, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) - a * pars[2])
}

#' Flat prior for GP parameters (\eqn{log \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 2.
#'   GP parameters (\eqn{\sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gp_flat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]))
}

#' Flat prior for GP parameters (\eqn{\sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 2.
#'   GP parameters (\eqn{\sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @param upper A positive numeric scalar.  The upper endpoint of the GP
#'   distribution.
#' @return The log of the prior density.
#' @export
gp_flatflat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0,
                        upper = NULL) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  if (!is.null(upper)) {
    if (pars[2] > min(0, max_xi) | -pars[1] / pars[2] > upper) {
      return(-Inf)
    }
  }
  return(0)
}

#' Jeffreys prior for GP parameters (\eqn{\sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 2.
#'   GP parameters (\eqn{\sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gp_jeffreys <- function(pars, min_xi = -1/2, max_xi = Inf, trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) - log(1 + pars[2]) - log(1 + 2 * pars[2]) / 2)
}

#' Beta-type prior for GP shape parameter \eqn{\xi}
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 2.
#'   GP parameters (\eqn{\sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param pq A numeric vector of length 2.
#'   See \code{\link{set_prior}} for details.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gp_beta <- function(pars, min_xi = -1 / 2, max_xi = 1 / 2, pq = c(6, 9),
                    trendsd = 0) {
  if (pars[1] <= 0 | pars[2] < min_xi | pars[2] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[1]) + (pq[1] - 1) * log(pars[2] - min_xi) +
           (pq[2] - 1) * log(max_xi - pars[2]))
}

# ================================= GEV priors =================================

gev_prior <- function(prior=c("norm", "loglognorm", "mdi", "flat", "flatflat",
                              "beta", "prob", "quant"), ...) {
  prior <- match.arg(prior)
  temp <- list(prior = paste("gev_", prior, sep=""), ...)
  # For v1.2.0 ...
  # Add default value for min_xi in the MDI prior if it isn't present.
  if (prior == "mdi" & is.null(temp$min_xi)) {
    temp$min_xi = -1
  }
  # Add default values for min_xi and max_xi in the beta prior if they
  # aren't present.
  if (prior == "beta" & is.null(temp$min_xi)) {
    temp$min_xi = -1/2
  }
  if (prior == "beta" & is.null(temp$max_xi)) {
    temp$max_xi = 1/2
  }
  # Force min_xi and max_xi to be present in the list temp.
  temp$min_xi <- max(temp$min_xi, -Inf)
  temp$max_xi <- min(temp$max_xi, Inf)
  # ... add default value for a in MDI prior if it isn't present.
  if (prior == "mdi" & is.null(temp$a)) {
    temp$a <- 0.5772156649015323
  }
  # ... add default value for pq in beta-type prior if it isn't present.
  if (prior == "beta" & is.null(temp$pq)) {
    temp$pq = c(6, 9)
  }
  # Check for unused hyperparameter names and drop them
  hpar_vec <- switch(prior, norm = c("mean", "cov"),
                     loglognorm = c("mean", "cov"), mdi = "a", flat = NULL,
                     beta = "pq", prob = c("quant", "alpha"),
                     quant = c("prob", "shape", "scale"))
  hpar_vec <- c(hpar_vec, "min_xi", "max_xi")
  temp <- hpar_drop(temp, hpar_vec)
  # Check for problems with min_xi and/or max_xi
  if (!is.null(temp$min_xi) & !is.null(temp$max_xi)) {
    if (temp$min_xi >= temp$max_xi)
        stop("min_xi must be less than max_xi")
  }
  if (!is.null(temp$min_xi)) {
    if (prior == "mdi" & is.infinite(temp$min_xi)) {
      stop("If min_xi=-Inf then the MDI posterior is not proper: min_xi must
            be finite.")
    }
  }
  # Check admissibility of hyperparameters
  if (prior == "norm" | prior == "loglognorm") {
    mean <- temp$mean
    cov <- temp$cov
    if (length(mean) != 3 | mode(mean) != "numeric")
        stop("mean must be a numeric vector of length three")
    if (is.null(cov))
        stop("cov must be supplied")
    if (!is.matrix(cov) | any(dim(cov) != 3) | mode(cov) != "numeric")
        stop("cov must be a symmetric three by three matrix")
    if (any(abs(cov - t(cov)) > .Machine$double.eps ^ 0.5))
        warning("cov may not be symmetric")
    eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
    if (any(eg <= 0))
        warning("cov may not be positive definite")
    icov <- solve(cov)
    icov <- icov[row(icov) >= col(icov)]
    temp$cov <- NULL
    temp <- c(temp, list(icov=icov))
  }
  if (prior == "mdi" & !is.null(temp$a)) {
    a <- temp$a
    if (length(a) != 1 | !is.numeric(a) | a < 0)
        stop("a must be a non-negative numeric vector of length 1")
  }
  if (prior == "beta" & !is.null(temp$pq)) {
    pq <- temp$pq
    if (length(pq) != 2 | !is.numeric(pq) | any(pq <= 0) )
        stop("pq must be a non-negative numeric vector of length 2")
  }
  if (prior == "prob") {
    if (is.null(temp$quant)) {
      stop("quant must be supplied when prior = prob")
    }
    if (is.null(temp$alpha)) {
      stop("alpha must be supplied when prior = prob")
    }
    if (length(temp$quant) != 3 | mode(temp$quant) != "numeric") {
      stop("quant must be a numeric vector of length three")
    }
    if (length(temp$alpha) != 4 | mode(temp$alpha) != "numeric") {
      stop("alpha must be a numeric vector of length four")
    }
    # Sort quant to ensure that the quantiles are increasing.
    temp$quant <- sort(temp$quant, decreasing = FALSE)
  }
  if (prior == "quant") {
    if (is.null(temp$prob)) {
      stop("prob must be supplied when prior = quant")
    }
    if (is.null(temp$shape)) {
      stop("shape must be supplied when prior = quant")
    }
    if (is.null(temp$scale)) {
      stop("scale must be supplied when prior = quant")
    }
    if (length(temp$prob) != 3 | mode(temp$prob) != "numeric") {
      stop("prob must be a numeric vector of length three")
    }
    if (length(temp$shape) != 3 | mode(temp$shape) != "numeric") {
      stop("shape must be a numeric vector of length three")
    }
    if (length(temp$scale) != 3 | mode(temp$scale) != "numeric") {
      stop("scale must be a numeric vector of length three")
    }
    # Sort prob to ensure that the corresponding quantiles are increasing.
    temp$prob <- sort(temp$prob, decreasing = TRUE)
  }
  # Add trendsd to the prior list so that the prior will work with the
  # evdbayes function posterior().
  temp$trendsd <- 0
  return(temp)
}

# ------------------------------ specific GEV priors ------------------------- #

#' Trivariate normal prior for GEV parameters (\eqn{\mu, log \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param mean A numeric vector of length 3.  Prior mean.
#' @param icov A 3x3 numeric matrix.
#'   The inverse of the prior covariance matrix.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_norm <- function(pars, mean, icov, min_xi = -Inf, max_xi = Inf,
                     trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  pars[2] <- log(pars[2])
  cpar <- pars - mean
  ld <- icov[1] * cpar[1] ^ 2 + 2 * icov[2] * cpar[1] * cpar[2] +
    2 * icov[3] * cpar[1] * cpar[3] + icov[4] * cpar[2] ^ 2 +
    2 * icov[5] * cpar[2] * cpar[3] + icov[6] * cpar[3] ^ 2
  return(-ld / 2 - pars[2])
}

#' Trivariate normal prior for GEV parameters (\eqn{log \mu, log \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param mean A numeric vector of length 3.  Prior mean.
#' @param icov A 3x3 numeric matrix.
#'   The inverse of the prior covariance matrix.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_loglognorm <- function(pars, mean, icov, min_xi = -Inf, max_xi = Inf,
                           trendsd = 0) {
  if (pars[1] <= 0 | pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  pars[1] <- log(pars[1])
  pars[2] <- log(pars[2])
  cpar <- pars - mean
  ld <- icov[1] * cpar[1] ^ 2 + 2 * icov[2] * cpar[1] * cpar[2] +
    2 * icov[3] * cpar[1] * cpar[3] + icov[4] * cpar[2] ^ 2 +
    2 * icov[5] * cpar[2] * cpar[3] + icov[6] * cpar[3] ^ 2
  return(-ld / 2 - pars[2] - pars[1])
}

#' Maximal data information (MDI) prior for GEV parameters
#' (\eqn{\mu, \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param a A numeric scalar.  The default value, Euler's constant, gives the
#'   MDI prior.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_mdi <- function(pars, a=0.5772156649015323, min_xi=-1, max_xi=Inf,
                    trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[2]) - a * pars[3])
}

#' Flat prior for GEV parameters (\eqn{\mu, log \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_flat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[2]))
}

#' Flat prior for GEV parameters (\eqn{\mu, \sigma, \xi})
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#'   Must not be \code{-Inf} because this results in an improper posterior.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_flatflat <- function(pars, min_xi = -Inf, max_xi = Inf, trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(0)
}

#' Beta-type prior for GEV shape parameter \eqn{\xi}
#'
#' For information about this and other priors see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param pq A numeric vector of length 2.
#'   See \code{\link{set_prior}} for details.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @export
gev_beta <- function(pars, min_xi = -1 / 2, max_xi = 1 / 2, pq = c(6, 9),
                     trendsd = 0) {
  if (pars[2] <= 0 | pars[3] < min_xi | pars[3] > max_xi) {
    return(-Inf)
  }
  return(-log(pars[2]) + (pq[1] - 1) * log(pars[3] - min_xi) +
           (pq[2] - 1) * log(max_xi - pars[3]))
}

# ================================ hpar_drop ===============================

hpar_drop <- function(x_list, hpar_vec) {
  # Check for unused hyperparameter names and drop them
  #
  # Args:
  #   x_list   : A list. A list to define a prior of class "evprior".
  #   hpar_vec : A character vector.  The names of prior hyperparameters.
  #
  # Returns:
  #   The input list x_list with any unused hyperparameters removed.
  #
  to_drop <- 1 + which(!(names(x_list[-1]) %in% hpar_vec))
  if (length(to_drop) == 1) {
    warning("This user-supplied argument is unused and has been dropped:",
            immediate. = TRUE)
  }
  if (length(to_drop) > 1) {
    warning("The following user-supplied arguments are unused and have been
            dropped:", immediate. = TRUE)
  }
  if (length(to_drop) > 0) {
    cat(names(x_list[to_drop]), "\n")
    x_list <- x_list[-to_drop]
  }
  return(x_list)
}

# ============================== set_bin_prior ================================

#' Construction of a prior distribution for a binomial probability \eqn{p}
#'
#' Constructs a prior distribution for use as the argument \code{bin_prior} in
#' \code{\link{rpost}} or in \code{\link{binpost}}.  The user can choose
#' from a list of in-built priors or specify their own prior function,
#' returning the \strong{log} of the prior density, using an R function
#' and arguments for hyperparameters.
#'
#' @param prior Either
#' \itemize{
#'   \item {An R function that returns the value of the log of the prior
#'   density (see \strong{Examples}), or}
#'   \item {A character string giving the name of the prior for \eqn{p}.
#'     See \strong{Details} for a list of priors available.}
#' }
#' @param ... Further arguments to be passed to the user-supplied or in-built
#'   prior function.  For the latter this is only relevant if
#'   \code{prior = "beta"}, when \code{ab} can be passed. See \strong{Details}.
#' @details
#'   \strong{Binomial priors.} The names of the binomial priors set using
#'   \code{bin_prior} are:
#' \itemize{
#'   \item{\code{"jeffreys"}: the \emph{Jeffreys} beta(1/2, 1/2) prior.}
#'   \item{\code{"laplace"}: the \emph{Bayes-Laplace} beta(1, 1) prior.}
#'   \item{\code{"haldane"}: the \emph{Haldane} beta(0, 0) prior.}
#'   \item{\code{"beta"}: a beta(\eqn{\alpha, \beta}) prior.  The argument
#'     \code{ab} is a vector containing \code{c}(\eqn{\alpha, \beta}).
#'     The default is \code{ab = c(1, 1)}.}
#'   \item{\code{"mdi"}: the MDI prior
#'     \eqn{\pi(p) = 1.6186 p^p (1-p)^{1-p}}{\pi(p) = 1.6186 p^p (1-p)^(1-p)},
#'         for \eqn{0 < p < 1.}}
#'   \item{\code{"northrop"}: the improper prior
#'     \eqn{\pi(p)=\{-\ln(1-p)\}^{-1}(1-p)^{-1}}{\pi(p)=1 / [ -ln(1-p) (1-p) ]},
#'         for \eqn{0 < p < 1.}}
#' }
#' Apart from the last two priors these are all beta distributions.
#' @return A list of class \code{"binprior"}.  The first component is the
#'   name of the input prior.  Apart from the MDI prior this will be "beta",
#'   in which case the other component of the list is a vector of length two
#'   giving the corresponding values of the beta parameters.
#' @seealso \code{\link{binpost}} for sampling from a binomial posterior
#'   distribution.
#' @examples
#' bp <- set_bin_prior(prior = "jeffreys")
#'
#' # Setting the Jeffreys prior by hand
#' beta_prior_fn <- function(p, ab) {
#'   return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
#' }
#' jeffreys <- set_bin_prior(beta_prior_fn, ab = c(1 / 2, 1 / 2))
#' @export
set_bin_prior <- function(prior = c("jeffreys", "laplace", "haldane", "beta",
                                    "mdi", "northrop"), ...) {
  # If prior is a function then return it in the required format.
  if (is.function(prior)) {
    temp <- list(prior = prior, ...)
    return(structure(temp, class = "binprior"))
  }
  if (class(prior) == "externalptr") {
    stop("A user-supplied prior must be specified using an R function")
  }
  # Otherwise, call the appropriate function to set the prior with name prior.
  prior <- match.arg(prior)
  temp <- list(prior = paste("bin_", prior, sep=""), ...)
  # Check for unused hyperparameter names and drop them
  if (prior == "beta") {
    hpar_vec <- "ab"
  } else {
    hpar_vec <- NULL
  }
  temp <- hpar_drop(temp, hpar_vec)
  # Check admissibility of hyperparameters
  if (prior == "beta" & !is.null(temp$ab)) {
    ab <- temp$ab
    if (length(ab) != 2 | !is.numeric(ab) | any(ab <= 0) )
      stop("ab must be a non-negative numeric vector of length 2")
  }
  if (prior == "beta" & is.null(temp$ab)) {
    temp$ab <- c(1, 1)
  }
  # Jeffreys, Laplace and Haldane prior are all beta distributions.
  temp$ab <- switch(prior, jeffreys = c(1/2, 1/2), laplace = c(1, 1),
               haldane = c(0, 0), beta = temp$ab)
  if (prior %in% c("jeffreys", "laplace", "haldane")) {
    temp$prior <- paste("bin_", "beta", sep = "")
  }
  return(structure(temp, class = "binprior"))
}

