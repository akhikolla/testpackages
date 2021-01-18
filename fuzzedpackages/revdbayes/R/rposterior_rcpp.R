# =========================== rpost_rcpp ===========================
#
#' Random sampling from extreme value posterior distributions
#'
#' Uses the \code{\link[rust]{ru_rcpp}} function in the
#' \code{\link[rust]{rust}} package to simulate from the posterior distribution
#' of an extreme value model.
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param model A character string.  Specifies the extreme value model.
#' @param data  Sample data, of a format appropriate to the value of
#'   \code{model}.
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
#' @param ... Further arguments to be passed to \code{\link[rust]{ru_rcpp}}.
#'   Most importantly \code{trans} and \code{rotate} (see \strong{Details}),
#'   and perhaps \code{r}, \code{ep}, \code{a_algor}, \code{b_algor},
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
#'   \code{\link[rust]{ru_rcpp}}.  In addition this list contains
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
#'     \item{\code{logf_rho_args}:} A list of arguments to the (transformed)
#'       target log-density.
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
#' @seealso \code{\link{rpost}} for posterior simulation without using
#'   the Rcpp package.
#' @seealso \code{\link{plot.evpost}}, \code{\link{summary.evpost}} and
#'   \code{\link{predict.evpost}} for the S3 \code{plot}, \code{summary}
#'   and \code{predict} methods for objects of class \code{evpost}.
#' @seealso \code{\link[rust]{ru_rcpp}} in the \code{\link[rust]{rust}}
#'   package for details of the arguments that can be passed to
#'   \code{ru_rcpp} and the form of the object returned by \code{rpost_rcpp}.
#' @seealso \code{\link[rust]{find_lambda}} in the
#'   \code{\link[rust]{rust}} package is used inside \code{rpost} to set the
#'   Box-Cox transformation parameter lambda when the \code{trans = "BC"}
#'   argument is given.
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
#' # GP model
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost_rcpp(n = 1000, model = "gp", prior = fp, thresh = u,
#'                   data = gom)
#' plot(gpg)
#'
#' # GP model, user-defined prior (same prior as the previous example)
#' ptr_gp_flat <- create_prior_xptr("gp_flat")
#' p_user <- set_prior(prior = ptr_gp_flat, model = "gp", min_xi = -1)
#' gpg <- rpost_rcpp(n = 1000, model = "gp", prior = p_user, thresh = u,
#'                   data = gom)
#' plot(gpg)
#'
#' # Binomial-GP model
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' bp <- set_bin_prior(prior = "jeffreys")
#' bgpg <- rpost_rcpp(n = 1000, model = "bingp", prior = fp, thresh = u,
#'                    data = gom, bin_prior = bp)
#' plot(bgpg, pu_only = TRUE)
#' plot(bgpg, add_pu = TRUE)
#'
#' # Setting the same binomial (Jeffreys) prior by hand
#' beta_prior_fn <- function(p, ab) {
#'   return(stats::dbeta(p, shape1 = ab[1], shape2 = ab[2], log = TRUE))
#' }
#' jeffreys <- set_bin_prior(beta_prior_fn, ab = c(1 / 2, 1 / 2))
#' bgpg <- rpost_rcpp(n = 1000, model = "bingp", prior = fp, thresh = u,
#'                    data = gom, bin_prior = jeffreys)
#' plot(bgpg, pu_only = TRUE)
#' plot(bgpg, add_pu = TRUE)
#'
#' # GEV model
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0, 0, 0), cov = mat)
#' gevp  <- rpost_rcpp(n = 1000, model = "gev", prior = pn, data = portpirie)
#' plot(gevp)
#'
#' # GEV model, user-defined prior (same prior as the previous example)
#' mat <- diag(c(10000, 10000, 100))
#' ptr_gev_norm <- create_prior_xptr("gev_norm")
#' pn_u <- set_prior(prior = ptr_gev_norm, model = "gev", mean = c(0, 0, 0),
#'                   icov = solve(mat))
#' gevu <- rpost_rcpp(n = 1000, model = "gev", prior = pn_u, data = portpirie)
#' plot(gevu)
#'
#' # GEV model, informative prior constructed on the probability scale
#' pip  <- set_prior(quant = c(85, 88, 95), alpha = c(4, 2.5, 2.25, 0.25),
#'   model = "gev", prior = "prob")
#' ox_post <- rpost_rcpp(n = 1000, model = "gev", prior = pip, data = oxford)
#' plot(ox_post)
#'
#' # PP model
#' pf <- set_prior(prior = "flat", model = "gev", min_xi = -1)
#' ppr <- rpost_rcpp(n = 1000, model = "pp", prior = pf, data = rainfall,
#'                   thresh = 40, noy = 54)
#' plot(ppr)
#'
#' # PP model, user-defined prior (same prior as the previous example)
#' ptr_gev_flat <- create_prior_xptr("gev_flat")
#' pf_u <- set_prior(prior = ptr_gev_flat, model = "gev", min_xi = -1,
#'                   max_xi = Inf)
#' ppru <- rpost_rcpp(n = 1000, model = "pp", prior = pf_u, data = rainfall,
#'                    thresh = 40, noy = 54)
#' plot(ppru)
#'
#' # PP model, informative prior constructed on the quantile scale
#' piq <- set_prior(prob = 10^-(1:3), shape = c(38.9, 7.1, 47),
#'                  scale = c(1.5, 6.3, 2.6), model = "gev", prior = "quant")
#' rn_post <- rpost_rcpp(n = 1000, model = "pp", prior = piq, data = rainfall,
#'                       thresh = 40, noy = 54)
#' plot(rn_post)
#'
#' # OS model
#' mat <- diag(c(10000, 10000, 100))
#' pv <- set_prior(prior = "norm", model = "gev", mean = c(0, 0, 0), cov = mat)
#' osv <- rpost_rcpp(n = 1000, model = "os", prior = pv, data = venice)
#' plot(osv)
#' @export
rpost_rcpp <- function(n, model = c("gev", "gp", "bingp", "pp", "os"), data,
                       prior, ..., nrep = NULL, thresh = NULL, noy = NULL,
                       use_noy = TRUE, npy = NULL, ros= NULL,
                       bin_prior = structure(list(prior = "bin_beta",
                                             ab = c(1 / 2, 1 / 2),
                                             class = "binprior")),
                       init_ests = NULL, mult = 2, use_phi_map = FALSE) {
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
                     use_noy = use_noy, ros = ros)
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
    ds$n_raw <- NULL
    temp_bin <- binpost(n = n, prior = bin_prior, ds_bin = ds_bin)
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
  # v1.2.0
  #
  # Create a pointer to the log-posterior function.
  #
  if (class(prior$prior) == "externalptr") {
    prior_type <- "user"
    post_ptr <- switch(model,
                       gp = gp_logpost_xptr("gp_user"),
                       gev = gev_logpost_xptr("gev_user"),
                       pp = pp_logpost_xptr("gev_user"),
                       os = os_logpost_xptr("gev_user")
    )
  } else {
    prior_type <- prior$prior
    post_ptr <- switch(model,
                       gp = gp_logpost_xptr(prior$prior),
                       gev = gev_logpost_xptr(prior$prior),
                       pp = pp_logpost_xptr(prior$prior),
                       os = os_logpost_xptr(prior$prior)
    )
  }
  #
  # Combine lists ds (data) and prior (details of prior) into one list.
  #
  for_post <- c(ds, prior)
  # For the OS model add the largest sample size over all order statistics.
  # (only used by the C++ function cpp_os_loglik)
  if (model == "os") {
    for_post <- c(for_post, list(nmax = n_check))
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
  init_check <- calc_init_logpost(model = model, prior_type = prior_type,
                                  init = init, for_post = for_post)
  #
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
    init_check <- calc_init_logpost(model = model, prior_type = prior_type,
                                    init = init_ests, for_post = for_post)
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
    for_ru_rcpp <- c(list(logf = post_ptr, pars = for_post), fr,
                     list(init = init, n = n), ru_args)
    temp <- do.call(rust::ru_rcpp, for_ru_rcpp)
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
  # v1.2.0
  # Create a pointers to the C++ phi_to_theta function and the C++ function
  # to evaluate the log-posterior after transformation from theta to phi.
  #
  phi_to_theta_ptr <- phi_to_theta_xptr(model)
  cpp_logpost_phi <- set_logpost_phi(model = model, prior_type = prior_type)
  #
  # Set which_lam: indices of the parameter vector that are Box-Cox transformed.
  #
  which_lam <- set_which_lam(model = model)
  #
  if (use_phi_map) {
    temp <- stats::optim(init_phi, cpp_logpost_phi,
                  control = list(parscale = se_phi, fnscale = -1),
                  pars = for_post)
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
  user_args <- switch(model,
    gp = list(xm = ds$xm),
    gev = list(x1 = ds$x1, xm = ds$xm),
    os = list(x1 = ds$x1, xm = ds$xm),
    pp = list(thresh = ds$thresh, xm = ds$xm)
  )
  for_find_lambda <- c(list(logf = post_ptr, pars = for_post), min_max_phi,
                       list(d = fr$d, which_lam = which_lam),
                       find_lambda_args, list(phi_to_theta = phi_to_theta_ptr,
                       user_args = user_args))
  min_max_phi$min_phi[3] <- 0.1
  lambda <- do.call(rust::find_lambda_rcpp, for_find_lambda)
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
  for_ru_rcpp <- c(list(logf = post_ptr, pars = for_post, lambda = lambda),
                   fr, list(n = n), ru_args)
  temp <- do.call(rust::ru_rcpp, for_ru_rcpp)
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

calc_init_logpost <- function(model, prior_type, init, for_post) {
  if (prior_type == "user") {
    init_check <- switch(
      model,
      gp = gp_user_logpost(x = init, pars = for_post),
      gev = gev_user_logpost(x = init, pars = for_post),
      os = os_user_logpost(x = init, pars = for_post),
      pp = pp_user_logpost(x = init, pars = for_post)
    )
    return(init_check)
  }
  if (model == "gp") {
    init_check <- switch(
      prior_type,
      gp_mdi = gp_mdi_logpost(x = init, pars = for_post),
      gp_norm = gp_norm_logpost(x = init, pars = for_post),
      gp_flat = gp_flat_logpost(x = init, pars = for_post),
      gp_flatflat = gp_flatflat_logpost(x = init, pars = for_post),
      gp_jeffreys = gp_jeffreys_logpost(x = init, pars = for_post),
      gp_beta = gp_beta_logpost(x = init, pars = for_post))
  }
  if (model == "gev") {
    init_check <- switch(
      prior_type,
      gev_mdi = gev_mdi_logpost(x = init, pars = for_post),
      gev_norm = gev_norm_logpost(x = init, pars = for_post),
      gev_loglognorm = gev_loglognorm_logpost(x = init, pars = for_post),
      gev_flat = gev_flat_logpost(x = init, pars = for_post),
      gev_flatflat = gev_flatflat_logpost(x = init, pars = for_post),
      gev_beta = gev_beta_logpost(x = init, pars = for_post),
      gev_prob = gev_prob_logpost(x = init, pars = for_post),
      gev_quant = gev_quant_logpost(x = init, pars = for_post))
  }
  if (model == "os") {
    init_check <- switch(
      prior_type,
      gev_mdi = os_mdi_logpost(x = init, pars = for_post),
      gev_norm = os_norm_logpost(x = init, pars = for_post),
      gev_loglognorm = os_loglognorm_logpost(x = init, pars = for_post),
      gev_flat = os_flat_logpost(x = init, pars = for_post),
      gev_flatflat = os_flatflat_logpost(x = init, pars = for_post),
      gev_beta = os_beta_logpost(x = init, pars = for_post),
      gev_prob = os_prob_logpost(x = init, pars = for_post),
      gev_quant = os_quant_logpost(x = init, pars = for_post))
  }
  if (model == "pp") {
    init_check <- switch(
      prior_type,
      gev_mdi = pp_mdi_logpost(x = init, pars = for_post),
      gev_norm = pp_norm_logpost(x = init, pars = for_post),
      gev_loglognorm = pp_loglognorm_logpost(x = init, pars = for_post),
      gev_flat = pp_flat_logpost(x = init, pars = for_post),
      gev_flatflat = pp_flatflat_logpost(x = init, pars = for_post),
      gev_beta = pp_beta_logpost(x = init, pars = for_post),
      gev_prob = pp_prob_logpost(x = init, pars = for_post),
      gev_quant = pp_quant_logpost(x = init, pars = for_post))
  }
  return(init_check)
}

set_logpost_phi <- function(model, prior_type) {
  if (model == "gp") {
    cpp_logpost_phi <- switch(prior_type,
                              gp_mdi = gp_mdi_logpost_phi,
                              gp_norm = gp_norm_logpost_phi,
                              gp_flat = gp_flat_logpost_phi,
                              gp_flatflat = gp_flatflat_logpost_phi,
                              gp_jeffreys = gp_jeffreys_logpost_phi,
                              gp_beta = gp_beta_logpost_phi,
                              user = gp_user_logpost_phi)
  } else if (model == "gev") {
    cpp_logpost_phi <- switch(prior_type,
                              gev_mdi = gev_mdi_logpost_phi,
                              gev_norm = gev_norm_logpost_phi,
                              gev_loglognorm = gev_loglognorm_logpost_phi,
                              gev_flat = gev_flat_logpost_phi,
                              gev_flatflat = gev_flatflat_logpost_phi,
                              gev_beta = gev_beta_logpost_phi,
                              gev_prob = gev_prob_logpost_phi,
                              gev_quant = gev_quant_logpost_phi,
                              user = gev_user_logpost_phi)
  } else if (model == "os") {
    cpp_logpost_phi <- switch(prior_type,
                              gev_mdi = os_mdi_logpost_phi,
                              gev_norm = os_norm_logpost_phi,
                              gev_loglognorm = os_loglognorm_logpost_phi,
                              gev_flat = os_flat_logpost_phi,
                              gev_flatflat = os_flatflat_logpost_phi,
                              gev_beta = os_beta_logpost_phi,
                              gev_prob = os_prob_logpost_phi,
                              gev_quant = os_quant_logpost_phi,
                              user = os_user_logpost_phi)
  } else if (model == "pp") {
    cpp_logpost_phi <- switch(prior_type,
                              gev_mdi = pp_mdi_logpost_phi,
                              gev_norm = pp_norm_logpost_phi,
                              gev_loglognorm = pp_loglognorm_logpost_phi,
                              gev_flat = pp_flat_logpost_phi,
                              gev_flatflat = pp_flatflat_logpost_phi,
                              gev_beta = pp_beta_logpost_phi,
                              gev_prob = pp_prob_logpost_phi,
                              gev_quant = pp_quant_logpost_phi,
                              user = pp_user_logpost_phi)
  }
  return(cpp_logpost_phi)
}
