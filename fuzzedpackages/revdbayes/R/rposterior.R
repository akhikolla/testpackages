# ================================ rpost ==================================== #
#
#' Random sampling from extreme value posterior distributions
#'
#' Uses the \code{\link[rust]{ru}} function in the \code{\link[rust]{rust}}
#' package to simulate from the posterior distribution of an extreme value
#' model.
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param model A character string.  Specifies the extreme value model.
#' @param data  Sample data, of a format appropriate to the value of
#'   \code{model}..
#'   \itemize{
#'     \item {"gp"} {A numeric vector of threshold excesses or raw data.}
#'     \item {"bingp"} {A numeric vector of raw data.}
#'     \item {"gev"} {A numeric vector of block maxima.}
#'     \item {"pp"} {A numeric vector of raw data.}
#'     \item {"os"} {A numeric matrix or data frame. Each row should contain
#'       the largest order statistics for a block of data.  These need not
#'       be ordered: they are sorted inside \code{rpost}. If a block
#'       contains fewer than \code{dim(as.matrix(data))[2]} order statistics
#'       then the corresponding row should be padded by \code{NA}s. If
#'       \code{ros} is supplied then only the largest \code{ros} values in
#'       each row are used.  If a vector is supplied then this is converted
#'       to a matrix with one column.  This is equivalent to using
#'       \code{model = "gev"}.}
#'   }
#' @param prior A list specifying the prior for the parameters of the extreme
#'   value model, created by \code{\link{set_prior}}.
#' @param nrep A numeric scalar.  If \code{nrep} is not \code{NULL} then
#'   \code{nrep} gives the number of replications of the original dataset
#'   simulated from the posterior predictive distribution.
#'   Each replication is based on one of the samples from the posterior
#'   distribution.  Therefore, \code{nrep} must not be greater than \code{n}.
#'   In that event \code{nrep} is set equal to \code{n}.
#'   Currently only implemented if \code{model = "gev"} or \code{"gp"} or
#'   \code{"bingp"} or \code{"pp"}, i.e. \emph{not} implemented if
#'   \code{model = "os"}.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#'   Only relevant when \code{model = "gp"}, \code{"pp"} or \code{"bingp"}.
#'   Must be supplied when \code{model = "pp"} or \code{"bingp"}.
#'   If \code{model = "gp"} and \code{thresh} is not supplied then
#'   \code{thresh = 0} is used and \code{data} should contain threshold
#'   excesses.
#' @param noy A numeric scalar. The number of blocks of observations,
#'   excluding any missing values.  A block is often a year.
#'   Only relevant, and must be supplied, if \code{model = "pp"}.
#' @param use_noy A logical scalar.  Only relevant if model is "pp".
#'   If \code{use_noy = FALSE} then sampling is based on a likelihood in
#'   which the number of blocks (years) is set equal to the number of threshold
#'   excesses, to reduce posterior dependence between the parameters
#'   (\href{https://doi.org/10.1214/10-AOAS333}{Wadsworth \emph{et al}. (2010)}).
#'   The sampled values are transformed back to the required parameterisation
#'   before returning them to the user.  If \code{use_noy = TRUE} then the
#'   user's value of \code{noy} is used in the likelihood.
#' @param npy A numeric scalar. The mean number of observations per year
#'   of data, after excluding any missing values, i.e. the number of
#'   non-missing observations divided by total number of years of non-missing
#'   data.
#'
#'   The value of \code{npy} does not affect any calculation in
#'   \code{rpost}, it only affects subsequent extreme value inferences using
#'   \code{predict.evpost}.  However, setting \code{npy} in the call to
#'   \code{rpost} avoids the need to supply \code{npy} when calling
#'   \code{predict.evpost}.  This is likely to be useful only when
#'   \code{model = bingp}. See the documentation of
#'   \code{\link{predict.evpost}} for further details.
#' @param ros A numeric scalar.  Only relevant when \code{model = "os"}. The
#'   largest \code{ros} values in each row of the matrix \code{data} are used
#'   in the analysis.
#' @param bin_prior A list specifying the prior for a binomial probability
#'   \eqn{p}, created by \code{\link{set_bin_prior}}.  Only relevant if
#'   \code{model = "bingp"}.  If this is not supplied then the Jeffreys
#'   beta(1/2, 1/2) prior is used.
#' @param bin_param A character scalar.  The argument \code{param} passed to
#'   \code{\link{binpost}}.  Only relevant if a user-supplied prior function
#'   is set using \code{\link{set_bin_prior}}.
#' @param init_ests A numeric vector.  Initial parameter estimates for search
#'   for the mode of the posterior distribution.
#' @param mult A numeric scalar.  The grid of values used to choose the Box-Cox
#'   transformation parameter lambda is based on the maximum a posteriori (MAP)
#'   estimate +/- mult x estimated posterior standard deviation.
#' @param use_phi_map A logical scalar. If trans = "BC" then \code{use_phi_map}
#'   determines whether the grid of values for phi used to set lambda is
#'   centred on the maximum a posterior (MAP) estimate of phi
#'   (\code{use_phi_map = TRUE}), or on the initial estimate of phi
#'   (\code{use_phi_map = FALSE}).
#' @param weights An optional numeric vector of weights by which to multiply
#'   the observations when constructing the log-likelihood.
#'   Currently only implemented for \code{model = "gp"} or
#'   \code{model = "bingp"}.
#'   In the latter case \code{bin_prior$prior} must be \code{"bin_beta"}.
#'   \code{weights} must have the same length as \code{data}.
#' @param ... Further arguments to be passed to \code{\link[rust]{ru}}.  Most
#'   importantly \code{trans} and \code{rotate} (see \strong{Details}), and
#'   perhaps \code{r}, \code{ep}, \code{a_algor}, \code{b_algor},
#'   \code{a_method}, \code{b_method}, \code{a_control}, \code{b_control}.
#'   May also be used to pass the arguments \code{n_grid} and/or \code{ep_bc}
#'   to \code{\link[rust]{find_lambda}}.
#' @details
#' \emph{Generalised Pareto (GP)}: \code{model = "gp"}.  A model for threshold
#'   excesses.  Required arguments: \code{n}, \code{data} and \code{prior}.
#'   If \code{thresh} is supplied then only the values in \code{data} that
#'   exceed \code{thresh} are used and the GP distribution is fitted to the
#'   amounts by which those values exceed \code{thresh}.
#'   If \code{thresh} is not supplied then the GP distribution is fitted to
#'   all values in \code{data}, in effect \code{thresh = 0}.
#'   See also \code{\link{gp}}.
#'
#' \emph{Binomial-GP}: \code{model = "bingp"}.  The GP model for threshold
#'   excesses supplemented by a binomial(\code{length(data)}, \eqn{p})
#'   model for the number of threshold excesses.  See
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   for details.  Currently, the GP and binomial parameters are assumed to
#'   be independent \emph{a priori}.
#'
#' \emph{Generalised extreme value (GEV) model}: \code{model = "gev"}.  A
#'   model for block maxima.  Required arguments: \code{n}, \code{data},
#'   \code{prior}.  See also \code{\link{gev}}.
#'
#' \emph{Point process (PP) model}: \code{model = "pp"}. A model for
#'   occurrences of threshold exceedances and threshold excesses.  Required
#'   arguments: \code{n}, \code{data}, \code{prior}, \code{thresh} and
#'   \code{noy}.
#'
#' \emph{r-largest order statistics (OS) model}: \code{model = "os"}.
#'   A model for the largest order statistics within blocks of data.
#'   Required arguments: \code{n}, \code{data}, \code{prior}.  All the values
#'   in \code{data} are used unless \code{ros} is supplied.
#'
#' \emph{Parameter transformation}.  The scalar logical arguments (to the
#'   function \code{\link{ru}}) \code{trans} and \code{rotate} determine,
#'   respectively, whether or not Box-Cox transformation is used to reduce
#'   asymmetry in the posterior distribution and rotation of parameter
#'   axes is used to reduce posterior parameter dependence.  The default
#'   is \code{trans = "none"} and \code{rotate = TRUE}.
#'
#'   See the \href{https://CRAN.R-project.org/package=revdbayes}{Introducing revdbayes vignette}
#'   for further details and examples.
#'
#' @return An object (list) of class \code{"evpost"}, which has the same
#'   structure as an object of class \code{"ru"} returned from
#'   \code{\link[rust]{ru}}.
#'   In addition this list contains
#'   \itemize{
#'     \item{\code{model}:} The argument \code{model} to \code{rpost}
#'       detailed above.
#'     \item{\code{data}:} The content depends on \code{model}:
#'       if \code{model = "gev"} then this is the argument \code{data} to
#'       \code{rpost} detailed above, with missing values removed;
#'       if \code{model = "gp"} then only the values that lie above the
#'       threshold are included; if \code{model = "bingp"} or
#'       \code{model = "pp"} then the input data are returned
#'       but any value lying below the threshold is set to \code{thresh};
#'       if \code{model = "os"} then the order statistics used are returned
#'       as a single vector.
#'     \item{\code{prior}:} The argument \code{prior} to \code{rpost}
#'       detailed above.
#'   }
#'   If \code{nrep} is not \code{NULL} then this list also contains
#'   \code{data_rep}, a numerical matrix with \code{nrep} rows.  Each
#'   row contains a replication of the original data \code{data}
#'   simulated from the posterior predictive distribution.
#'   If \code{model = "bingp"} or \code{"pp"} then the rate of threshold
#'   exceedance is part of the inference.  Therefore, the number of values in
#'   \code{data_rep} that lie above the threshold varies between
#'   predictive replications (different rows of \code{data_rep}).
#'   Values below the threshold are left-censored at the threshold, i.e. they
#'   are set at the threshold.
#'
#'   If \code{model == "pp"} then this list also contains the argument
#'     \code{noy} to \code{rpost} detailed above.
#'   If the argument \code{npy} was supplied then this list also contains
#'   \code{npy}.
#'
#'   If \code{model == "gp"} or \code{model == "bingp"} then this list also
#'     contains the argument \code{thresh} to \code{rpost} detailed above.
#'
#'   If \code{model == "bingp"} then this list also contains
#'   \itemize{
#'     \item{\code{bin_sim_vals}:} {An \code{n} by 1 numeric matrix of values
#'       simulated from the posterior for the binomial
#'       probability \eqn{p}}
#'     \item{\code{bin_logf}:} {A function returning the log-posterior for
#'       \eqn{p}.}
#'     \item{\code{bin_logf_args}:} {A list of arguments to \code{bin_logf}.}
#'   }
#' @seealso \code{\link{set_prior}} for setting a prior distribution.
#' @seealso \code{\link{rpost_rcpp}} for faster posterior simulation using
#'   the Rcpp package.
#' @seealso \code{\link{plot.evpost}}, \code{\link{summary.evpost}} and
#'   \code{\link{predict.evpost}} for the S3 \code{plot}, \code{summary}
#'   and \code{predict} methods for objects of class \code{evpost}.
#' @seealso \code{\link[rust]{ru}} and \code{\link[rust]{ru_rcpp}} in the
#'   \code{\link[rust]{rust}} package for details of the arguments that can
#'   be passed to \code{ru} and the form of the object returned by
#'   \code{rpost}.
#' @seealso \code{\link[rust]{find_lambda}} and
#'   \code{\link[rust]{find_lambda_rcpp}} in the \code{\link[rust]{rust}}
#'   package is used inside \code{rpost} to set the Box-Cox transformation
#'   parameter lambda when the \code{trans = "BC"} argument is given.
#' @seealso \code{\link[evdbayes]{posterior}} for sampling from an extreme
#'   value posterior using the evdbayes package.
#' @references Coles, S. G. and Powell, E. A. (1996) Bayesian methods in
#'   extreme value modelling: a review and new developments.
#'   \emph{Int. Statist. Rev.}, \strong{64}, 119-136.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
#' @references Stephenson, A. (2016) Bayesian Inference for Extreme Value
#'   Modelling. In \emph{Extreme Value Modeling and Risk Analysis: Methods and
#'   Applications}, edited by D. K. Dey and J. Yan, 257-80. London:
#'   Chapman and Hall. \url{https://doi.org/10.1201/b19721}
#'   value posterior using the evdbayes package.
#' @references Wadsworth, J. L., Tawn, J. A. and Jonathan, P. (2010)
#'   Accounting for choice of measurement scale in extreme value modeling.
#'  \emph{The Annals of Applied Statistics}, \strong{4}(3), 1558-1578.
#'   \url{https://doi.org/10.1214/10-AOAS333}
#' @examples
#' \donttest{
#' # GP model
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost(n = 1000, model = "gp", prior = fp, thresh = u, data = gom)
#' plot(gpg)
#'
#' # Binomial-GP model
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' bp <- set_bin_prior(prior = "jeffreys")
#' bgpg <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u, data = gom,
#'               bin_prior = bp)
#' plot(bgpg, pu_only = TRUE)
#' plot(bgpg, add_pu = TRUE)
#'
#' # Setting the same binomial (Jeffreys) prior by hand
#' beta_prior_fn <- function(p, ab) {
#'   return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
#' }
#' jeffreys <- set_bin_prior(beta_prior_fn, ab = c(1 / 2, 1 / 2))
#' bgpg <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u, data = gom,
#'               bin_prior = jeffreys)
#' plot(bgpg, pu_only = TRUE)
#' plot(bgpg, add_pu = TRUE)
#'
#' # GEV model
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0, 0, 0), cov = mat)
#' gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)
#' plot(gevp)
#'
#' # GEV model, informative prior constructed on the probability scale
#' pip  <- set_prior(quant = c(85, 88, 95), alpha = c(4, 2.5, 2.25, 0.25),
#'                   model = "gev", prior = "prob")
#' ox_post <- rpost(n = 1000, model = "gev", prior = pip, data = oxford)
#' plot(ox_post)
#'
#' # PP model
#' pf <- set_prior(prior = "flat", model = "gev", min_xi = -1)
#' ppr <- rpost(n = 1000, model = "pp", prior = pf, data = rainfall,
#'              thresh = 40, noy = 54)
#' plot(ppr)
#'
#' # PP model, informative prior constructed on the quantile scale
#' piq <- set_prior(prob = 10^-(1:3), shape = c(38.9, 7.1, 47),
#'                  scale = c(1.5, 6.3, 2.6), model = "gev", prior = "quant")
#' rn_post <- rpost(n = 1000, model = "pp", prior = piq, data = rainfall,
#'                  thresh = 40, noy = 54)
#' plot(rn_post)
#'
#' # OS model
#' mat <- diag(c(10000, 10000, 100))
#' pv <- set_prior(prior = "norm", model = "gev", mean = c(0, 0, 0), cov = mat)
#' osv <- rpost(n = 1000, model = "os", prior = pv, data = venice)
#' plot(osv)
#' }
#' @export
rpost <- function(n, model = c("gev", "gp", "bingp", "pp", "os"), data, prior,
                  ..., nrep = NULL, thresh = NULL, noy = NULL, use_noy = TRUE,
                  npy = NULL, ros= NULL,
                  bin_prior = structure(list(prior = "bin_beta",
                                             ab = c(1 / 2, 1 / 2),
                                             class = "binprior")),
                  bin_param = "logit",
                  init_ests = NULL, mult = 2, use_phi_map = FALSE,
                  weights = NULL) {
  #
  model <- match.arg(model)
  save_model <- model
  # Check that the prior is compatible with the model.
  # If an evdbayes prior has been set make the "model" attribute of the
  # prior equal to "gev"
  if (is.character(prior$prior)) {
    if (prior$prior %in% c("dprior.norm", "dprior.loglognorm", "dprior.quant",
                         "dprior.prob")) {
      attr(prior, "model") <- "gev"
    }
  }
  prior_model <- attr(prior, "model")
  if (prior_model %in% c("pp", "os")) {
    prior_model <- "gev"
  }
  check_model <- model
  if (model %in% c("gev", "pp", "os")) {
    check_model <- "gev"
  }
  if (model %in% c("gp", "bingp")) {
    check_model <- "gp"
  }
  if (prior_model != check_model) {
    stop("The prior is not compatible with the model.")
  }
  # ---------- Create list of additional arguments to the likelihood --------- #
  # Check that any required arguments to the likelihood are present.
  if (model == "gp" & is.null(thresh)) {
    thresh <- 0
    warning("threshold thresh was not supplied so thresh = 0 is used",
            immediate. = TRUE)
  }
  if (model == "bingp" & is.null(thresh)) {
    stop("threshold thresh must be supplied when model is bingp")
  }
  if (model == "pp" & is.null(thresh)) {
    stop("threshold thresh must be supplied when model is pp")
  }
  if (model == "pp" & is.null(noy)) {
    stop("number of years noy must be supplied when model is pp")
  }
  # ------------------------------ Process the data ------------------------- #
  # Remove missings, extract sample summaries, shift and scale the data to have
  # approximate location 0 and scale 1. This avoids numerical problems that may
  # result if the posterior scales of the parameters are very different.
  #
  ds <- process_data(model = model, data = data, thresh = thresh, noy = noy,
                     use_noy = use_noy, ros = ros, weights = weights)
  #
  # If model = "bingp" then extract sufficient statistics for the binomial
  # model, and remove n_raw from ds because it isn't used in the GP
  # log-likelihood.  Sample from the posterior distribution for the binomial
  # probability p.  Then set model = "gp" because we have dealt with the "bin"
  # bit of "bingp".
  #
  add_binomial <- FALSE
  if (model == "bingp") {
    ds_bin <- list()
    ds_bin$m <- ds$m
    ds_bin$n_raw <- ds$n_raw
    if (is.null(weights)) {
      temp_bin <- binpost(n = n, prior = bin_prior, ds_bin = ds_bin,
                          param = bin_param)
    } else {
      ds_bin$sf <- ds$sf
      ds_bin$w <- ds$binw
      ds$sf <- ds$binw <- NULL
      temp_bin <- wbinpost(n = n, prior = bin_prior, ds_bin = ds_bin)
    }
    ds$n_raw <- NULL
    add_binomial <- TRUE
    model <- "gp"
  }
  #
  n_check <- ds$n_check
  if (n_check < 10) {
    warning(paste(
      "With a small sample size of", n_check,
      "it may be that optimisations will fail."), immediate. = TRUE,
            noBreaks. = TRUE)
  }
  ds$n_check <- NULL
  # Check that if one of the in-built improper priors is used then the sample
  # size is sufficiently large to produce a proper posterior distribution.
  if (is.character(prior$prior)) {
    check_sample_size(prior_name = prior$prior, n_check = n_check)
  }
  # ----------- Extract min_xi and max_xi from prior (if supplied) ---------- #
  #
  min_xi <- ifelse(is.null(prior$min_xi), -Inf, prior$min_xi)
  max_xi <- ifelse(is.null(prior$max_xi), +Inf, prior$max_xi)
  #
  # -------------------------- Set up log-posterior ------------------------- #
  #
  if (model == "gp") {
    logpost <- function(pars, ds) {
      if (is.null(weights)) {
        loglik <- do.call(gp_loglik, c(list(pars = pars), ds))
      } else {
        loglik <- do.call(gp_wloglik, c(list(pars = pars), ds))
      }
      if (is.infinite(loglik)) return(loglik)
      logprior <- do.call(prior$prior, c(list(pars), prior[-1]))
      return(loglik + logprior)
    }
  }
  if (model == "gev") {
    logpost <- function(pars, ds) {
      loglik <- do.call(gev_loglik, c(list(pars = pars), ds))
      if (is.infinite(loglik)) return(loglik)
      logprior <- do.call(prior$prior, c(list(pars), prior[-1]))
      return(loglik + logprior)
    }
  }
  if (model == "os") {
    logpost <- function(pars, ds) {
      loglik <- do.call(os_loglik, c(list(pars = pars), ds))
      if (is.infinite(loglik)) return(loglik)
      logprior <- do.call(prior$prior, c(list(pars), prior[-1]))
      return(loglik + logprior)
    }
  }
  if (model == "pp") {
    if (ds$noy == noy) {
      logpost <- function(pars, ds) {
        loglik <- do.call(pp_loglik, c(list(pars = pars), ds))
        if (is.infinite(loglik)) return(loglik)
        logprior <- do.call(prior$prior, c(list(pars), prior[-1]))
        return(loglik + logprior)
      }
    } else {
      log_rat <- log(ds$noy / noy)
      logpost <- function(pars, ds) {
        loglik <- do.call(pp_loglik, c(list(pars = pars), ds))
        if (is.infinite(loglik)) return(loglik)
        noy_pars <- change_pp_pars(pars, in_noy = ds$noy, out_noy = noy)
        logprior <- do.call(prior$prior, c(list(noy_pars), prior[-1]))
        return(loglik + logprior + pars[3] * log_rat)
      }
    }
  }
  #
  # ---------------- set some admissible initial estimates ------------------ #
  #
  temp <- switch(model,
                 gp = do.call(gp_init, ds),
                 gev = do.call(gev_init, ds),
                 pp = do.call(pp_init, ds),
                 os = do.call(os_init, ds)
  )
  init <- temp$init
  #
  # If init is not admissible set xi = 0 and try again
  #
  init_check <- logpost(pars = init, ds = ds)
  if (is.infinite(init_check)) {
    temp <- switch(model,
                   gp = do.call(gp_init, c(ds, xi_eq_zero = TRUE)),
                   gev = do.call(gev_init, c(ds, xi_eq_zero = TRUE)),
                   pp = do.call(pp_init, c(ds, xi_eq_zero = TRUE)),
                   os = do.call(os_init, c(ds, xi_eq_zero = TRUE))
    )
    init <- temp$init
  }
  se <- temp$se
  init_phi <- temp$init_phi
  se_phi <- temp$se_phi
  #
  # Use user's initial estimates (if supplied and if admissible)
  #
  if (!is.null(init_ests)) {
    # If model = "pp" and use_noy = FALSE, so that we are working with a
    # parameterisation in which the number of blocks is equal to the number
    # of threshold excesses, then transform the user-supplied initial
    # estimates from the number of blocks = noy parameterisation to the number
    # of blocks = number of threshold excesses.
    if (model == "pp" & use_noy == FALSE) {
      init_ests <- change_pp_pars(init_ests, in_noy = noy, out_noy = ds$n_exc)
    }
    init_check <- logpost(par = init_ests, ds = ds)
    if (!is.infinite(init_check)) {
      init <- init_ests
      init_phi <- switch(model,
                         gp = do.call(gp_init, c(ds, list(init_ests
                                                          = init_ests))),
                         gev = do.call(gev_init, c(ds, list(init_ests
                                                            = init_ests))),
                         pp = do.call(pp_init, c(ds, list(init_ests
                                                            = init_ests))),
                         os = do.call(os_init, c(ds, list(init_ests
                                                          = init_ests)))
      )
    }
  }
  #
  # Extract arguments to be passed to ru()
  #
  ru_args <- list(...)
  #
  # Set default values for trans and rotate if they have not been supplied.
  #
  if (is.null(ru_args$trans)) ru_args$trans <- "none"
  if (is.null(ru_args$rotate)) ru_args$rotate <- TRUE
  #
  # Create list of objects to send to function ru()
  #
  fr <- create_ru_list(model = model, trans = ru_args$trans,
                       rotate = ru_args$rotate, min_xi = min_xi,
                       max_xi = max_xi)
  #
  # In anticipation of the possibility of inclusion of a one parameter model
  # in the future.
  #
  if (fr$d == 1) {
    ru_args$rotate <- FALSE
  }
  #
  if (ru_args$trans == "none") {
    # Only set a_control$parscale if it hasn't been supplied and if a_algor
    # will be "optim" in ru()
    if (is.null(ru_args$a_control$parscale)) {
      if (is.null(ru_args$a_algor) & fr$d > 1) {
        ru_args$a_control$parscale <- se
      }
      if (!is.null(ru_args$a_algor)) {
        if (ru_args$a_algor == "optim") {
          ru_args$a_control$parscale <- se
        }
      }
    }
    #
    # Set ru_args$n_grid and ru_args$ep_bc to NULL just in case they have been
    # specified in ...
    #
    ru_args$n_grid <- NULL
    ru_args$ep_bc <- NULL
    for_ru <- c(list(logf = logpost, ds = ds), fr, list(init = init, n = n),
                ru_args)
    temp <- do.call(rust::ru, for_ru)
    #
    # If model == "pp" and the sampling parameterisation is not equal to that
    # required by the user then transform to the required parameterisation.
    #
    if (model == "pp") {
      if (ds$noy != noy) {
        temp$sim_vals <- t(apply(temp$sim_vals, 1, FUN = change_pp_pars,
                                  in_noy = ds$noy, out_noy = noy))
      }
    }
    # If model was "bingp" then add the binomial posterior simulated values.
    if (add_binomial) {
      temp$bin_sim_vals <- matrix(temp_bin$bin_sim_vals, ncol = 1)
      colnames(temp$bin_sim_vals) <- "p[u]"
      temp$bin_logf <- temp_bin$bin_logf
      temp$bin_logf_args <- temp_bin$bin_logf_args
    }
    class(temp) <- "evpost"
    temp$model <- save_model
    if (save_model == "gp") {
      temp$data <- ds$data + thresh
    } else if (save_model == "bingp") {
      temp$data <- ds$data + thresh
      temp$data <- c(temp$data, rep(thresh, ds_bin$n_raw - length(temp$data)))
    } else if (save_model == "pp") {
      temp$data <- ds$data
      temp$data <- c(temp$data, rep(thresh, ds$m - length(temp$data)))
    } else {
      temp$data <- ds$data
    }
    if (save_model == "pp") {
      temp$noy <- noy
    }
    if (save_model %in% c("gp", "bingp")) {
      temp$thresh <- thresh
    }
    temp$npy <- npy
    temp$prior <- prior
    if (!is.null(nrep)) {
      if (save_model == "os") {
        warning("model = ``os'' so nrep has been ignored.")
      }
      if (nrep > n & save_model != "os") {
        nrep <- n
        warning("nrep has been set equal to n.")
      }
      if (save_model == "gev") {
        wr <- 1:nrep
        temp$data_rep <- replicate(ds$m, rgev(nrep, loc = temp$sim_vals[wr, 1],
                                              scale = temp$sim_vals[wr, 2],
                                              shape = temp$sim_vals[wr, 3]))
      }
      if (save_model == "gp") {
        wr <- 1:nrep
        temp$data_rep <- replicate(ds$m, rgp(nrep, loc = thresh,
                                            scale = temp$sim_vals[wr, 1],
                                            shape = temp$sim_vals[wr, 2]))
      }
      if (save_model == "bingp") {
        wr <- 1:nrep
        temp$data_rep <- matrix(thresh, nrow = nrep, ncol = ds_bin$n_raw)
        for (i in wr) {
          n_above <- stats::rbinom(1, ds_bin$n_raw, temp$bin_sim_vals[i])
          if (n_above > 0) {
            temp$data_rep[i, 1:n_above] <- rgp(n = n_above, loc = thresh,
                                               scale = temp$sim_vals[i, 1],
                                               shape = temp$sim_vals[i, 2])
          }
        }
      }
      if (save_model == "pp") {
        wr <- 1:nrep
        temp$data_rep <- matrix(thresh, nrow = nrep, ncol = length(temp$data))
        loc <- temp$sim_vals[wr, 1]
        scale <- temp$sim_vals[wr, 2]
        shape <- temp$sim_vals[wr, 3]
        mod_scale <- scale + shape * (thresh - loc)
        p_u <- noy * pu_pp(q = thresh, loc = loc, scale = scale,
                           shape = shape, lower_tail = FALSE) / ds$m
        for (i in wr) {
          n_above <- stats::rbinom(1, ds$m, p_u)
          if (n_above > 0) {
            temp$data_rep[i, 1:n_above] <- rgp(n = n_above, loc = thresh,
                                               scale = mod_scale[i],
                                               shape = shape[i])
          }
        }
      }
    }
    return(temp)
  }
  #
  # ----------------- If Box-Cox transformation IS required ----------------- #
  #
  # -------------------------- Define phi_to_theta -------------------------- #
  #
  if (model == "gp") {
    phi_to_theta <- function(phi) {
      c(phi[1], phi[2] - phi[1] / ds$xm)
    }
    log_j <- function(x) 0
  }
  if (model == "gev" | model == "os") {
    sr <- sqrt(ds$xm - ds$x1)
    phi_to_theta <- function(phi) {
      mu <- phi[1]
      xi <- (phi[3] - phi[2]) / sr
      sigma <- ((ds$xm - phi[1]) * phi[2] + (phi[1] - ds$x1) * phi[3]) / sr
      c(mu, sigma, xi)
    }
    log_j <- function(x) 0
  }
  if (model == "pp") {
    sr <- sqrt(ds$xm - ds$thresh)
    phi_to_theta <- function(phi) {
      mu <- phi[1]
      xi <- (phi[3] - phi[2]) / sr
      sigma <- ((ds$xm - phi[1]) * phi[2] + (phi[1] - ds$thresh) * phi[3]) / sr
      c(mu, sigma, xi)
    }
    log_j <- function(x) 0
  }
  #
  # Set which_lam: indices of the parameter vector that are Box-Cox transformed.
  #
  which_lam <- set_which_lam(model = model)
  #
  if (use_phi_map) {
    logpost_phi <- function(phi, ...) {
      logpost(par = phi_to_theta(phi), ...)
    }
    temp <- stats::optim(init_phi, logpost_phi,
                  control = list(parscale = se_phi, fnscale = -1), ds = ds)
    phi_mid <- temp$par
  } else {
    phi_mid <- init_phi
  }
  #
  min_max_phi <- set_range_phi(model = model, phi_mid = phi_mid,
                               se_phi = se_phi, mult = mult)
  #
  # Look for find_lambda arguments n_grid and ep_bc in ...
  #
  temp_args <- list(...)
  find_lambda_args <- list()
  if (!is.null(temp_args$n_grid)) find_lambda_args$n_grid <- temp_args$n_grid
  if (!is.null(temp_args$ep_bc)) find_lambda_args$n_grid <- temp_args$ep_bc
  # Find Box-Cox parameter lambda and initial estimate for psi.
  for_find_lambda <- c(list(logf = logpost, ds = ds), min_max_phi,
                       list(d = fr$d, which_lam = which_lam),
                       find_lambda_args,
                       list(phi_to_theta = phi_to_theta, log_j = log_j))
  lambda <- do.call(rust::find_lambda, for_find_lambda)
  #
  # Only set a_control$parscale if it hasn't been supplied and if a_algor
  # will be "optim" in ru()
  #
  if (is.null(ru_args$a_control$parscale)) {
    if (is.null(ru_args$a_algor) & fr$d > 1) {
      ru_args$a_control$parscale <- lambda$sd_psi
    }
    if (!is.null(ru_args$a_algor)) {
      if (ru_args$a_algor == "optim") {
        ru_args$a_control$parscale <- lambda$sd_psi
      }
    }
  }
  for_ru <- c(list(logf = logpost, ds = ds, lambda = lambda), fr,
              list(n = n, phi_to_theta = phi_to_theta, log_j = log_j),
              ru_args)
  temp <- do.call(rust::ru, for_ru)
  #
  # If model == "pp" and the sampling parameterisation is not equal to that
  # required by the user then transform to the required parameterisation.
  #
  if (model == "pp") {
    if (ds$noy != noy) {
      temp$sim_vals <- t(apply(temp$sim_vals, 1, FUN = change_pp_pars,
                               in_noy = ds$noy, out_noy = noy))
    }
  }
  # If model was "bingp" then add the binomial posterior simulated values.
  if (add_binomial) {
    temp$bin_sim_vals <- matrix(temp_bin$bin_sim_vals, ncol = 1)
    colnames(temp$bin_sim_vals) <- "p[u]"
    temp$bin_logf <- temp_bin$bin_logf
    temp$bin_logf_args <- temp_bin$bin_logf_args
  }
  class(temp) <- "evpost"
  temp$model <- save_model
  if (save_model == "gp") {
    temp$data <- ds$data + thresh
  } else if (save_model == "bingp") {
    temp$data <- ds$data + thresh
    temp$data <- c(temp$data, rep(thresh, ds_bin$n_raw - length(temp$data)))
  } else if (save_model == "pp") {
    temp$data <- ds$data
    temp$data <- c(temp$data, rep(thresh, ds$m - length(temp$data)))
  } else {
    temp$data <- ds$data
  }
  if (save_model == "pp") {
    temp$noy <- noy
  }
  if (save_model %in% c("gp", "bingp")) {
    temp$thresh <- thresh
  }
  temp$npy <- npy
  temp$prior <- prior
  if (!is.null(nrep)) {
    if (nrep > n) {
      nrep <- n
      warning("nrep has been set equal to n.")
    }
    if (save_model == "gev") {
      wr <- 1:nrep
      temp$data_rep <- replicate(ds$m, rgev(nrep, loc = temp$sim_vals[wr, 1],
                                            scale = temp$sim_vals[wr, 2],
                                            shape = temp$sim_vals[wr, 3]))
    }
    if (save_model == "gp") {
      wr <- 1:nrep
      temp$data_rep <- replicate(ds$m, rgp(nrep, loc = thresh,
                                           scale = temp$sim_vals[wr, 1],
                                           shape = temp$sim_vals[wr, 2]))
    }
    if (save_model == "bingp") {
      wr <- 1:nrep
      temp$data_rep <- matrix(thresh, nrow = nrep, ncol = ds_bin$n_raw)
      for (i in wr) {
        n_above <- stats::rbinom(1, ds_bin$n_raw, temp$bin_sim_vals[i])
        temp$data_rep[i, 1:n_above] <- rgp(n = n_above, loc = thresh,
                                           scale = temp$sim_vals[i, 1],
                                           shape = temp$sim_vals[i, 2])
      }
    }
    if (save_model == "pp") {
      wr <- 1:nrep
      temp$data_rep <- matrix(thresh, nrow = nrep, ncol = length(temp$data))
      loc <- temp$sim_vals[, 1]
      scale <- temp$sim_vals[, 2]
      shape <- temp$sim_vals[, 3]
      mod_scale <- scale + shape * (thresh - loc)
      p_u <- noy * (mod_scale / scale) ^ (-1 / shape) / ds$m
      for (i in wr) {
        n_above <- stats::rbinom(1, ds$m, p_u[i])
        temp$data_rep[i, 1:n_above] <- rgp(n = n_above, loc = thresh,
                                           scale = mod_scale[i],
                                           shape = shape[i])
      }
    }
  }
  return(temp)
}

# ================================ pu_pp ==================================== #

pu_pp <- function (q, loc = 0, scale = 1, shape = 0, lower_tail = TRUE){
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  len_loc <- length(loc)
  len_scale <- length(scale)
  len_shape <- length(shape)
  check_len <- c(len_loc, len_scale, len_shape)
  if (length(unique(check_len[check_len > 1])) > 1) {
    stop("loc, scale and shape have incompatible lengths.")
  }
  max_len <- max(check_len)
  loc <- rep(loc, length.out = max_len)
  scale <- rep(scale, length.out = max_len)
  shape <- rep(shape, length.out = max_len)
  len_q <- length(q)
  if (max_len != 1 & len_q != 1) {
    stop("If length(q) > 1 then scale and shape must have length 1.")
  }
  q <- (q - loc) / scale
  nn <- length(q)
  qq <- 1 + shape * q
  # co is a condition to ensure that calculations are only performed in
  # instances where the density is positive.
  co <- qq > 0 | is.na(qq)
  q <- q[co]
  qq <- qq[co]
  p <- numeric(nn)
  p[!co] <- 1
  # If shape is close to zero then base calculation on an approximation that
  # is linear in shape.
  if (len_q == 1) {
    shape <- shape[co]
    p[co] <- ifelse(abs(shape) < 1e-6, 1 - exp(-q + shape * q ^ 2 / 2),
                    1 - pmax(1 + shape * q, 0) ^ (-1 / shape))
  } else {
    if(abs(shape) < 1e-6) {
      p[co] <- 1 - exp(-q + shape * q ^ 2 / 2)
    } else {
      p[co] <- 1 - pmax(1 + shape * q, 0) ^ (-1 / shape)
    }
  }
  if (!lower_tail) {
    p <- 1 - p
  }
  return(p)
}

# =============================== binpost =================================== #

#' Random sampling from a binomial posterior distribution
#'
#' Samples from the posterior distribution of the probability \eqn{p}
#' of a binomial distribution.
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param prior A function to evaluate the prior, created by
#'   \code{\link{set_bin_prior}}.
#' @param ds_bin A numeric list.  Sufficient statistics for inference
#'   about a binomial probability \eqn{p}.  Contains
#' \itemize{
#'   \item {\code{n_raw} : number of raw observations}
#'   \item {\code{m} : number of threshold exceedances.}
#' }
#' @param param A character scalar.  Only relevant if \code{prior$prior} is a
#'   (user-supplied) R function.  \code{param} specifies the parameterization
#'   of the posterior distribution that \code{\link[rust]{ru}} uses for
#'   sampling.
#'
#'   If \code{param = "p"} the original parameterization \eqn{p} is
#'   used.
#'
#'   If \code{param = "logit"} (the default) then \code{\link[rust]{ru}}
#'   samples from the posterior for the logit of \eqn{p}, before transforming
#'   back to the \eqn{p}-scale.
#'
#'   The latter tends to make the optimizations involved in the
#'   ratio-of-uniforms algorithm more stable and to increase the probability
#'   of acceptance, but at the expense of slower function evaluations.
#' @details If \code{prior$prior == "bin_beta"} then the posterior for \eqn{p}
#'   is a beta distribution so \code{\link[stats:Beta]{rbeta}} is used to
#'   sample from the posterior.
#'
#'   If \code{prior$prior == "bin_mdi"} then
#'   rejection sampling is used to sample from the posterior with an envelope
#'   function equal to the density of a
#'   beta(\code{ds$m} + 1, \code{ds$n_raw - ds$m} + 1) density.
#'
#'   If \code{prior$prior == "bin_northrop"} then
#'   rejection sampling is used to sample from the posterior with an envelope
#'   function equal to the posterior density that results from using a
#'   Haldane prior.
#'
#'   If \code{prior$prior} is a (user-supplied) R function then
#'   \code{\link[rust]{ru}} is used to sample from the posterior using the
#'   generalised ratio-of-uniforms method.
#' @return An object (list) of class \code{"binpost"} with components
#'   \itemize{
#'     \item{\code{bin_sim_vals}:} {An \code{n} by 1 numeric matrix of values
#'       simulated from the posterior for the binomial
#'       probability \eqn{p}}
#'     \item{\code{bin_logf}:} {A function returning the log-posterior for
#'       \eqn{p}.}
#'     \item{\code{bin_logf_args}:} {A list of arguments to \code{bin_logf}.}
#'   }
#'   If \code{prior$prior} is a (user-supplied) R function then this list
#'   also contains \code{ru_object} the object of class \code{"ru"}
#'   returned by \code{\link[rust]{ru}}.
#' @seealso \code{\link{set_bin_prior}} for setting a prior distribution
#'   for the binomial probability \eqn{p}.
#' @examples
#' data(gom)
#' u <- quantile(gom, probs = 0.65)
#' ds_bin <- list()
#' ds_bin$n_raw <- length(gom)
#' ds_bin$m <- sum(gom > u)
#' bp <- set_bin_prior(prior = "jeffreys")
#' temp <- binpost(n = 1000, prior = bp, ds_bin = ds_bin)
#' graphics::hist(temp$bin_sim_vals, prob = TRUE)
#'
#' # Setting a beta prior (Jeffreys in this case) by hand
#' beta_prior_fn <- function(p, ab) {
#'   return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
#' }
#' jeffreys <- set_bin_prior(beta_prior_fn, ab = c(1 / 2, 1 / 2))
#' temp <- binpost(n = 1000, prior = jeffreys, ds_bin = ds_bin)
#' @export
binpost <- function(n, prior, ds_bin, param = c("logit", "p")) {
  param <- match.arg(param)
  n_success <- ds_bin$m
  n_failure <- ds_bin$n_raw - ds_bin$m
  if (is.character(prior$prior) && prior$prior == "bin_beta") {
    shape1 <- n_success + prior$ab[1]
    shape2 <- n_failure + prior$ab[2]
    bin_sim_vals <- stats::rbeta(n = n, shape1 = shape1, shape2 = shape2)
    bin_logf_args <- list(shape1 = shape1, shape2 = shape2)
    bin_logf_beta <- function(x, shape1 , shape2) {
      stats::dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)
    }
    temp <- list(bin_sim_vals = bin_sim_vals, bin_logf = bin_logf_beta,
                 bin_logf_args = bin_logf_args)
  }
  if (is.character(prior$prior) && prior$prior == "bin_mdi") {
    bin_sim_vals <- r_pu_MDI(n = n, n_success = n_success,
                             n_failure = n_failure)
    beta_const <- beta(n_success + 1, n_failure + 1)
    log_beta_const <- log(beta_const)
    MDI_post <- function(pu) {
      pu ^ (pu + n_success) * (1 - pu) ^ (1 - pu + n_failure) / beta_const
    }
    log_MDI_const <- log(stats::integrate(MDI_post, 0, 1)$value)
    bin_logf_args <- list(n_success = n_success, n_failure = n_failure,
                          log_MDI_const = log_MDI_const,
                          log_beta_const = log_beta_const)
    bin_logf_mdi <- function(pu, n_success, n_failure, log_MDI_const,
                             log_beta_const) {
      if (pu > 1 || pu < 0) return(NA)
      (pu + n_success) * log(pu) + (1 - pu + n_failure) * log(1 - pu) -
        log_MDI_const - log_beta_const
    }
    temp <- list(bin_sim_vals = bin_sim_vals, bin_logf = bin_logf_mdi,
                 bin_logf_args = bin_logf_args)
  }
  if (is.character(prior$prior) && prior$prior == "bin_northrop") {
    bin_sim_vals <- r_pu_N(n = n, n_success = n_success, n_failure = n_failure)
    beta_const <- beta(n_success, n_failure)
    log_beta_const <- log(beta_const)
    N_post <- function(pu) {
      pu ^ n_success * (1 - pu) ^ (n_failure - 1) / (-log(1 - pu)) / beta_const
    }
    log_N_const <- log(stats::integrate(N_post, 0, 1)$value)
    bin_logf_args <- list(n_success = n_success, n_failure = n_failure,
                          log_N_const = log_N_const,
                          log_beta_const = log_beta_const)
    bin_logf_N <- function(pu, n_success, n_failure, log_N_const,
                           log_beta_const) {
      if (pu > 1 || pu < 0) return(NA)
      n_success * log(pu) + (n_failure - 1) * log(1 - pu) - log(-log(1 - pu)) -
        log_N_const - log_beta_const
    }
    temp <- list(bin_sim_vals = bin_sim_vals, bin_logf = bin_logf_N,
                 bin_logf_args = bin_logf_args)
  }
  if (is.function(prior$prior)) {
    binlogpost <- function(pu, ...) {
      binloglik <- stats::dbeta(pu, shape1 = n_success + 1,
                                shape2 = n_failure + 1, log = TRUE)
      binlogprior <- do.call(prior$prior, c(list(pu), prior[-1]))
      return(binloglik + binlogprior)
    }
    binpostfn <- function(pu, ...) {
      exp(binlogpost(pu, ...))
    }
    log_user_const <- log(stats::integrate(binpostfn, 0, 1)$value)
    bininit <- n_success / (n_success + n_failure)
    for_ru <- c(list(logf = binlogpost), list(n = n, d =  1, init = bininit))
    if (param == "logit") {
      for_ru$trans <- "user"
      # Transformation from phi (logit(p)) to theta (p)
      phi_to_theta <- function(phi) {
        return(exp(phi) / (1 + exp(phi)))
      }
      # Log-Jacobian of the transformation from theta to phi, i.e. based on the
      # derivatives of phi with respect to theta
      log_j <- function(theta) {
        return(-log(theta) - log(1 - theta))
      }
      for_ru$phi_to_theta <- phi_to_theta
      for_ru$log_j <- log_j
      for_ru$init <- log(bininit / (1 - bininit))
    }
    for_ru$var_names <- "p"
    temp2 <- do.call(rust::ru, for_ru)
    bin_logf_user <- function(pu, ...) {
      binlogpost(pu, ...) - log_user_const
    }
    temp <- list()
    temp$bin_sim_vals <- temp2$sim_vals
    temp$bin_logf <- bin_logf_user
    temp$bin_logf_args <- c(list(n_success = n_success, n_failure = n_failure),
                          prior[-1])
    temp$ru_object <- temp2
  }
  class(temp) <- "binpost"
  return(temp)
}

# ================================ r_pu_MDI ================================= #

r_pu_MDI <- function(n, n_success, n_failure) {
  #
  # Uses rejection sampling to sample from the posterior distribution of
  # a binomial probability p, based on the MDI prior for p.
  #
  # Args:
  #   n         : A numeric scalar.  The sample size required.
  #   n_success : A numeric scalar.  The number of successes.
  #   n_failure : A numeric scalar.  The number of failures.
  # Returns:
  #   A numeric vector containing the n sampled values.
  #
  n_acc <- 0
  p_sim <- NULL
  while (n_acc < n) {
    p <- stats::rbeta(1, n_success + 1, n_failure + 1)
    u <- stats::runif(1)
    if (u <= p ^ p * (1 - p) ^ (1 - p)) {
      n_acc <- n_acc+1
      p_sim[n_acc] <- p
    }
  }
  return(p_sim)
}

# ================================= r_pu_N ================================== #

r_pu_N <- function(n, n_success, n_failure, alpha = 0, beta = 0) {
  #
  # Uses rejection sampling to sample from the posterior distribution of
  # a binomial probability p, based on the prior pi(p) that is proportional to
  #             1 / [ -ln(1 - p) * (1 - p) ], for 0 < p < 1.
  # We use the posterior under a Beta(alpha, beta) prior as the envelope.
  # The posterior under this prior is Beta(n_success + alpha, n_failure + beta).
  # The default case, alpha = beta = 0, corresponds to using the posterior
  # under a Haldane prior as the envelope.  A slight improvement in the
  # probability of acceptance results from (alpha, beta) = (0.040842 0.005571).
  #
  # Args:
  #   n         : A numeric scalar.  The sample size required.
  #   n_success : A numeric scalar.  The number of successes.
  #   n_failure : A numeric scalar.  The number of failures.
  #   alpha     : A numeric scalar in [0, 1). Argument shape1 of rbeta()          alpha
  #   beta      : A numeric scalar in [0, 1). Argument shape1 of rbeta()          alpha
  # Returns:
  #   A numeric vector containing the n sampled values.
  #
  n_acc <- 0
  p_sim <- NULL
  while (n_acc < n) {
    p <- stats::rbeta(1, n_success + alpha, n_failure + beta)
    u <- stats::runif(1)
    if (u <= -p / log(1 - p)) {
      n_acc <- n_acc+1
      p_sim[n_acc] <- p
    }
  }
  return(p_sim)
}

# =============================== wbinpost ================================== #

#' Random sampling from a binomial posterior distribution, using weights
#'
#' Samples from the posterior distribution of the probability \eqn{p}
#' of a binomial distribution.  User-supplied weights are applied to each
#' observation when constructing the log-likelihood.
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param prior A function to evaluate the prior, created by
#'   \code{\link{set_bin_prior}}.
#'   \code{prior$prior} must be \code{"bin_beta"}.
#' @param ds_bin A numeric list.  Sufficient statistics for inference
#'   about the binomial probability \eqn{p}.  Contains
#' \itemize{
#'   \item {\code{sf} : a logical vector of success (\code{TRUE}) and failure
#'     (\code{FALSE}) indicators.}
#'   \item {\code{w} : a numeric vector of length \code{length(sf)} containing
#'     the values by which to multiply the observations when constructing the
#'     log-likelihood.}
#' }
#' @details For \code{prior$prior == "bin_beta"} the posterior for \eqn{p}
#'   is a beta distribution so \code{\link[stats:Beta]{rbeta}} is used to
#'   sample from the posterior.
#' @return An object (list) of class \code{"binpost"} with components
#'   \itemize{
#'     \item{\code{bin_sim_vals}:} {An \code{n} by 1 numeric matrix of values
#'       simulated from the posterior for the binomial
#'       probability \eqn{p}}
#'     \item{\code{bin_logf}:} {A function returning the log-posterior for
#'       \eqn{p}.}
#'     \item{\code{bin_logf_args}:} {A list of arguments to \code{bin_logf}.}
#'   }
#' @seealso \code{\link{set_bin_prior}} for setting a prior distribution
#'   for the binomial probability \eqn{p}.
#' @examples
#' u <- quantile(gom, probs = 0.65)
#' ds_bin <- list(sf = gom > u, w = rep(1, length(gom)))
#' bp <- set_bin_prior(prior = "jeffreys")
#' temp <- wbinpost(n = 1000, prior = bp, ds_bin = ds_bin)
#' graphics::hist(temp$bin_sim_vals, prob = TRUE)
#' @export
wbinpost <- function(n, prior, ds_bin) {
  n_success <- ds_bin$m
  n_failure <- ds_bin$n_raw - ds_bin$m
  if (is.character(prior$prior) && prior$prior == "bin_beta") {
    shape1 <- sum(ds_bin$w * ds_bin$sf) + prior$ab[1]
    shape2 <- sum(ds_bin$w * (1 - ds_bin$sf)) + prior$ab[2]
    bin_sim_vals <- stats::rbeta(n = n, shape1 = shape1, shape2 = shape2)
    bin_logf_args <- list(shape1 = shape1, shape2 = shape2)
    bin_logf_beta <- function(x, shape1 , shape2) {
      stats::dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE)
    }
    temp <- list(bin_sim_vals = bin_sim_vals, bin_logf = bin_logf_beta,
                 bin_logf_args = bin_logf_args)
  }
  class(temp) <- "binpost"
  return(temp)
}
