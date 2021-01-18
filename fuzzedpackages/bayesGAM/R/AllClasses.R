
## ------------- bayesGAM Classes ----------------------------------------

setClass("distribution",
         slots=c(name = "character",
                 distnum = "integer",
                 num_params = "integer",
                 param_names = "character",
                 param_values = "numeric"))

valid_distribution <- function(object) {
  if (length(object@name) != 1) {
    paste("Invalid distribution name", object@name, "Should be character of length one")
  }
  if (length(object@param_names) != length(object@param_values)) {
    paste("length of param_names must match length of param_values")
  } else {
    TRUE
  }
}

setValidity("distribution", valid_distribution)

# constructor for distributions
dtr <- distribution <- function(name, distnum, param_names, param_values) {
  new("distribution",
      name=name,
      distnum=distnum,
      param_names=param_names,
      param_values=param_values)
}

# -----------------normal distribution class---------------------------- #
setClass("normalDistribution",
         contains = "distribution")

#' Constructor function for Normal priors
#' 
#' Used to specify Normal priors for \code{bayesGAM} models
#' @param param_values Numeric vector of length 2 for the mean and standard deviation parameters
#' @details For the \code{beta} and \code{a} parameters, the distribution is assumed to be unconstrained. 
#' For \code{eps} and \code{lambda}, the priors are half-normal with a support of strictly positive numbers. 
#' @export
#' @references Stan Development Team. 2018. Stan Modeling Language Users Guide and Reference Manual, Version 2.18.0
#' @examples
#' require(stats); require(graphics)
#' normal(c(0, 10))
normal <- function(param_values) {
  new("normalDistribution",
      param_values=param_values)
}


# -----------------Student t distribution class---------------------------- #

setClass("tDistribution",
         contains = "distribution")

#' Constructor function for Student-t priors
#' 
#' Used to specify student-t priors for \code{bayesGAM} models
#' @param param_values Numeric vector of length 3 for the degrees of freedom, location, and scale parameter. 
#' @details For the \code{beta} and \code{a} parameters, the distribution is assumed to be unconstrained. 
#' For \code{eps} and \code{lambda}, the priors are half-normal with a support of strictly positive numbers. 
#' @export
#' @references Stan Development Team. 2018. Stan Modeling Language Users Guide and Reference Manual, Version 2.18.0
#' @examples
#' require(stats); require(graphics)
#' st(c(3,0,1))
st <- function(param_values) {
  new("tDistribution",
      param_values=param_values)
}

# -----------------prior class---------------------------- #

# must be a list of distributions
setClass("prior",
         slots = c(famnum = "integer",
                   mixed = "logical",
                   beta="list",
                   eps = "list",
                   lambda = "list",
                   a = "list",
                   random_intercept = "logical"))

# -----------------GLM family class---------------------------- #

setClass("glmFamily",
         slots=c(famnum = "integer",
                 linknum = "integer",
                 famname = "character",
                 linkname = "character",
                 mixed = "logical"))

setClass("glmModel",
         slots = c(X = "matrix",
                   Z = "matrix",
                   Zint = "matrix",
                   Znp = "matrix", 
                   # Zlst = "list",
                   # Zarray = "array",
                   max_col = "integer",
                   y = "matrix",
                   p = "integer",
                   r = "integer",
                   q = "integer",
                   n = "integer",
                   has_intercept = "logical",
                   zvars = "integer",
                   names_beta = "character",
                   names_u = "character",
                   names_y = "character",
                   prior = "prior",
                   knots = "list",
                   basis = "character",
                   npdegree = "integer",
                   npargs = "list",
                   npterms = "character",
                   sub_form = "formula", 
                   random_intercept = "logical",
                   multresponse = "logical", 
                   offset = "numeric", 
                   call = "call"),
                   # beta = "list",
                   # eps = "list",
                   # lambda = "list"),
                 contains = "glmFamily")


# ----------------- GLM model fit class ---------------------------- #
# --------------- post model fitting -----------------------------#

#' Contains results from \code{rstan} as well as the design matrices and other data for the model.
#' 
#' Returns object generated from model fit by \code{bayesGAM}
#' 
#' @aliases bayesGAMfit
#' @param object Object of type \code{bayesGAMfit} which stores the results from \code{rstan}, design matrices, and other data for the model.
#' @slot results Object of type \code{stanfit} returned from calling \code{rstan::sampling}
#' @slot model Object of custom type \code{glmModel} with the data and input parameters passed to \code{rstan}
#' @slot offset Optionally numeric offset for the generalized additive model
#' @slot spcontrol List of control parameters for \code{bayesGAMfit}
#' @slot mcmcres Matrix of MCMC results for all chains, if plot data is stored
#' @slot pdata Dataframe for default plot method, if plot data is stored.
#' @rdname bayesGAMfit
#' @export
setClass("bayesGAMfit",
         slots = c(results = "stanfit",
                   model = "glmModel",
                   offset = "numeric",
                   spcontrol = "list",
                   mcmcres = "matrix", 
                   pdata = "data.frame"))


# --------------- smooth plot object -----------------------------#


setClass("smoothPlotObject",
         slots = c(xvars = "character",
                   betanms = "character",
                   unms = "character",
                   betavals = "numeric",
                   uvals = "numeric",
                   xvars_static = "character",
                   xvars_npargs = "list",
                   xvars_np = "character",
                   xvars_bivariate = "list",
                   xvars_basis = "character",
                   npdegree = "integer",
                   has_intercept = "logical",
                   random_intercept = "logical",
                   zvars = "integer",
                   Xorig ="matrix",
                   knots= "list",
                   linkname = "character",
                   names_y = "character",
                   results = "stanfit",
                   multresponse = "logical", 
                   mcmcres = "matrix"))

# ----------------------- posterior predict object -------------------------#

setClass("posteriorPredictObject", 
         slots = c(model = "glmModel", 
                   pp = "list", 
                   thetasamp = "matrix"))


# ----------------------- predict object -------------------------#

setClass("predictPlotObject", 
         slots = c(model = "glmModel", 
                   pp = "list"))

