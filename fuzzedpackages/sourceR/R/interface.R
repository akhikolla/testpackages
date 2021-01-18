#####################################################
# Name: interface.R                                 #
# Author: Poppy Miller <p.miller@lancaster.ac.uk>   #
# Created: 2016-06-15                               #
# Copyright: Poppy Miller 2016                      #
# Purpose: Source attribution model interface       #
#####################################################

#' Builds a HaldDP source attribution model
#'
#' @docType class
#' @name HaldDP
#' @importFrom R6 R6Class
#' @importFrom grDevices col2rgb colorRampPalette
#' @importFrom stats median
#' @import dplyr
#' @import tensorA
#' @import assertthat
#' @export
#'
#' @param y a \code{\link{Y}} object containing case observations
#' @param x an \code{\link{X}} object containing source observations
#' @param k a \code{\link{Prev}} object containing source prevalences
#' @param priors \code{priors} list with elements named \code{a_r}, \code{a_alpha}, \code{a_theta} and \code{b_theta},
#'   corresponding to the prior parameters for the \code{r}, \code{alpha}, and base
#'   distribution for the DP parameters respectively.
#'
#'   \tabular{lllll}{
#'   \emph{Parameter} \tab \emph{Prior Distribution} \tab \emph{Prior Parameters}\cr
#'   \code{a_r} \tab Dirichlet(concentration) \tab A single positive number or an \code{\link{X}} \cr
#'   \tab \tab  object containing the prior values for each source,\cr
#'   \tab \tab time and type. If a single number is supplied,\cr
#'   \tab \tab it will be used for all times, sources and types. \cr
#'
#'   \code{a_alpha} \tab Dirichlet(concentration) \tab A single positive number or an \code{\link{Alpha}} \cr
#'   \tab \tab  object containing the prior values for each source,\cr
#'   \tab \tab time and location. If a single number is supplied,\cr
#'   \tab \tab it will be used for all times, sources and locations. \cr
#'
#'   Type effects base \tab Gamma(shape, rate) \tab Single number for each of the shape (a_theta) and \cr
#'   distribution parameters \tab \tab rate (b_theta) of the Gamma base distribution.\cr
#'   }
#'
#' @param inits initial values for the mcmc algorithm. This is an optional list
#'   that may contain any of the following items: \code{alpha},\code{q}, and \code{r}.
#'
#'   \tabular{lll}{
#'   \emph{Parameter} \tab \emph{Description} \cr
#'   \code{r}
#'   \tab An object of type \code{\link{X}} giving the initial values for $R$ matrix,\cr
#'   \tab If not specified defaults to the element-wise maximum likelihood\cr
#'   \tab estimates of \code{r} from the source matrix.\cr
#'   Source effects (\code{alpha})
#'   \tab An object of type \code{\link{Alpha}} specifying alpha value for each source/time/location.\cr
#'   \tab If not specified, default initial values\cr
#'   \tab for the source effects are drawn from the prior distribution. \cr
#'   Type effects (\code{q})
#'   \tab An object of type \code{\link{Q}} giving the initial clustering and values for \eqn{q}\cr
#'   \tab If not specified, defaults to a single group with a theta value calculated as \cr
#'   \tab \eqn{\theta = sum(y_itl) / sum_l=1^L(sum_t=1^T(sum_i=1^n(sum_j=1^m(alpha_jtl * r_ijt * k_jt))))}. \cr
#'   \tab i.e. \eqn{theta = sum(y_itl) / sum(lambda_ijtl / theta)}
#'   }
#'
#' @param a_q the Dirichlet Process concentration parameter.
#'
#' @format Object of \code{\link{R6Class}} with methods for creating a HaldDP model,
#' running the model, and accessing and plotting the results.
#' @section Description:
#' This function fits a non-parametric Poisson source attribution model for human cases of
#' disease. It supports multiple types, sources, times and locations. The number of
#' human cases for each type, time and location follow a Poisson likelihood.
#' @return Object of \code{\link{HaldDP}} with methods for creating a HaldDP model,
#' running the model, and accessing and plotting the results.
#'
#' @section HaldDP Object Methods:
#' \describe{
#'   \item{\code{mcmc_params(n_iter = 1000, burn_in = 0, thin = 1,
#'   n_r = ceiling(private$nTypes * 0.2), update_schema = c('q','alpha','r'))}}{when called, sets the mcmc
#'   parameters.
#'
#'   \code{n_iter} sets the number of iterations returned (after removing
#'   \code{burn_in} and thinning results by \code{thin} i.e. a total of
#'   (n_iter * thin) + burn_in iterations are run)
#'
#'   \code{n_r} is a positive
#'   integer that sets the total number of \code{r_{ijtl}} parameters to be updated
#'   at each time-location-source combination (the default is 20 percent updated
#'   per iteration)
#'
#'   \code{update_schema} a character vector containing the parameters to update
#'   (any of '\code{q}','\code{alpha}','\code{r}').
#'   }
#'
#'   \item{\code{update(n_iter, append = TRUE)}}{when called, updates the \code{HaldDP}
#'   model by running \code{n_iter} iterations.
#'
#'   If missing \code{n_iter}, the \code{n_iter} last set using \code{mcmc_params()}
#'   or \code{update()} is used.
#'
#'   \code{append}
#'   is a logical value which determines whether the next \code{n_iter} iterations
#'   are appended to any previous iterations, or overwrites them. When
#'   \code{append = TRUE}, the starting values are the last iteration and no
#'   \code{burn_in} is removed. Running the model for the first time, or changing any
#'   model or fitting parameters will set \code{append = FALSE}. }
#'
#'   \item{\code{get_data}}{returns a list containing the human data \code{y}
#'   (an array y[types, times, locations]), the source data \code{X} (an array X[types, sources, times]),
#'   the prevalence data (an array k[sources, times]), the type names, source names,
#'   time names, location names and number of different types, sources, times and locations.
#'   }
#'
#'   \item{\code{get_priors}}{returns a list containing the DP concentration
#'   parameter \code{a_q}, and the priors (R6 class with members named \code{a_alpha}
#'   (members are array \code{a_alpha[sources, times, locations]}), \code{a_r} (an array
#'   \code{a_r[types, sources, times]}), \code{a_theta} and \code{b_theta}).}
#'
#'   \item{\code{get_inits}}{returns an R6 class holding the initial values
#'   (members are \code{alpha} (an array \code{alpha[sources, times, locations]}),
#'   \code{theta} (an array \code{theta[types, iters]}), \code{s} (an array
#'   \code{s[types, iters]}), and \code{r} (an array \code{r[types, sources, times]})).}
#'
#'   \item{\code{get_mcmc_params}}{returns a list of fitting parameters (\code{n_iter},
#'   \code{append}, \code{burn_in}, \code{thin}, \code{update_schema} (R6 class with members
#'   \code{alpha}, \code{q}, \code{r})).}
#'
#'   \item{\code{get_acceptance}}{returns an R6 class containing the acceptance
#'   rates for each parameter (members are \code{alpha} (an array \code{alpha[sources, times, locations]}),
#'   and \code{r} (an array \code{r[types, sources, times]})).}
#'
#'   \item{\code{extract(params = c("alpha", "q", "s", "r", "lambda_i", "xi", "xi_prop"),
#'   times = NULL, locations = NULL, sources = NULL, types = NULL, iters = NULL,
#'   flatten = FALSE, drop = TRUE)}}{returns a list contining a subset of the parameters
#'   (determined by the \code{params} vector, \code{times}, \code{locations}, \code{sources}, \code{types} and \code{iters}).
#'
#'   If \code{flatten} is set to \code{TRUE}, it returns a dataframe with 1 column per
#'   parameter, otherwise it returns a list containing \code{params} containing a
#'   subset of the following arrays: \code{alpha[Sources, Times, Locations, iters]}, \code{q[Types, iters]},
#'   \code{s[Types, iters]}, \code{r[Types, Sources, Times, iters]},
#'   \code{lambda_i[Types, Times, Locations, iters]},
#'   \code{xi[Sources, Times, Locations, iters]}.
#'
#'   \code{drop}
#'   determines whether to delete the dimensions of an array which have only one
#'   level when \code{flatten = FALSE}.}
#'
#'   \item{\code{summary(alpha = 0.05, params = c("alpha", "q", "s", "r", "lambda_i",
#'   "xi" ,"xi_prop"), times = NULL, locations = NULL, sources = NULL,
#'   types = NULL, iters = NULL, flatten = FALSE, drop = TRUE, CI_type = "chen-shao")}}{
#'   returns a list contining the
#'   median and credible intervals for a subset of the parameters. The default credible
#'   interval type are Chen-Shao (\code{"chen-shao"}) highest posterior density intervals (alternatives
#'   are \code{"percentiles"} and \code{"spin"}).
#'   See \code{extract} for details on the subsetting. \code{xi_prop} returns the
#'   proportion of cases attributed to each source \code{j} and is calculated by dividing
#'   each iteration of \code{lambda_{jtl}} values by their sum within each time \code{t}
#'   and location \code{l}.}
#'
#'   \item{\code{plot_heatmap(iters, cols = c("blue","white"), hclust_method = "complete")}}{
#'   Creates a dendrogram and heatmap for the type effect groupings (\code{s} parameter
#'   in the model). This uses the heatmap.2 function from gplots.
#'
#'   \code{iters} is a vector containing the iterations to be used in constructing
#'   the graph. Default is all iterations in posterior.
#'
#'   \code{hclust_method} allows the user to select the method used by \code{stats::hclust} to
#'   cluster the type effect groupings \code{s}.
#'
#'   \code{cols} gives the colours for completely dissimilar (dissimilarity value
#'   of 1), and identical (dissimilarity value of 0). All other values will be in
#'   between the two chosen colours. See ?colorRampPalette for more details..}
#' }
#'
#' @section Details:
#' \describe{
#' This function fits a source attribution model for human cases of disease.
#' It supports multiple types, sources, times and locations. The number of human cases
#' for each type, time and location follows a Poisson or Negative Binomial likelihood.
#' \emph{Model}
#' \deqn{y_{itl}\sim\textsf{Poisson}(\lambda_{itl})}
#' where
#' \deqn{\lambda_{itl}=\sum_{j=1}^{m}\lambda_{ijtl}=q_{k(i)}\sum_{j=1}^{m}(r_{ijt}\cdot k_{j}\cdot alpha_{jtl})}
#'
#' The parameters are defined as follows:
#' \deqn{a_{jtl}} is the unknown source effect for source \eqn{j}, time \eqn{t}, location \eqn{l}
#' \deqn{q_{s(i)}} is the unknown type effect for type \eqn{i} in group \eqn{s}.
#' \deqn{x_{ij}} is the known number of positive samples for each source \eqn{j} type\eqn{i} combination
#' \deqn{n_{ij}} is the known total number of samples for each source \eqn{j} type \eqn{i} combination
#' \deqn{k_{j}} is the fixed prevalence in source (i.e. the number of positive samples
#' divided by the number of negative samples) \eqn{j}
#' \deqn{r_{ijt}}  is the unknown relative occurrence of type \eqn{i} on source \eqn{j}.
#'
#' \emph{Priors}
#' \deqn{r_{.jt}\sim Dirichlet(a\_r_{1jt},..., a\_r_{njt})}
#' \deqn{a_{tl}\sim Dirichlet(a\_alpha_{1tl},..., a\_alpha_{mtl})}
#' \deqn{q\sim DP(a_q, Gamma(a_{theta},b_{theta}))}
#' }
#'
#' @references Chen, M.-H. and Shao, Q.-M. (1998). Monte Carlo estimation of Bayesian
#' credible and HPD intervals, \emph{Journal of Computational and Graphical Statistics}, 7.
#' @references Liu Y, Gelman A, Zheng T (2015). Simulation-efficient shortest probability
#' intervals, \emph{Statistics and Computing}.
#' @author Chris Jewell and Poppy Miller \email{p.miller at lancaster.ac.uk}
#'
#' @examples
#'
#' #### Format data using Y, X, and Prev functions #############################
#' ## Input data must be in long format
#' y <- Y(                      # Cases
#'   data = sim_SA$cases,
#'   y = "Human",
#'   type = "Type",
#'   time = "Time",
#'   location = "Location"
#' )
#'
#' x <- X(                      # Sources
#'   data = sim_SA$sources,
#'   x = "Count",
#'   type = "Type",
#'   time = "Time",
#'   source = "Source"
#' )
#'
#' k <- Prev(                   # Prevalences
#'   data = sim_SA$prev,
#'   prev = "Value",
#'   time = "Time",
#'   source = "Source"
#' )
#'
#' #### Create Dirichlet(1) priors #############################################
#'
#' ## Create alpha prior data frame
#' prior_alpha_long <- expand.grid(
#'   Source   = unique(sim_SA$sources$Source),
#'   Time     = unique(sim_SA$sources$Time),
#'   Location = unique(sim_SA$cases$Location),
#'   Alpha    = 1
#' )
#' # Use the Alpha() constructor to specify alpha prior
#' prior_alpha <- Alpha(
#'   data     = prior_alpha_long,
#'   alpha    = 'Alpha',
#'   source   = 'Source',
#'   time     = 'Time',
#'   location = 'Location'
#' )
#'
#' ## Create r prior data frame
#' prior_r_long <- expand.grid(
#'   Type   = unique(sim_SA$sources$Type),
#'   Source = unique(sim_SA$sources$Source),
#'   Time   = unique(sim_SA$sources$Time),
#'   Value  = 0.1
#' )
#' # Use X() constructor to specify r prior
#' prior_r <- X(
#'   data   = prior_r_long,
#'   x      = 'Value',
#'   type   = 'Type',
#'   time   = 'Time',
#'   source = 'Source'
#' )
#'
#' ## Pack all priors into a list
#' priors <- list(
#'   a_theta = 0.01,
#'   b_theta = 0.00001,
#'   a_alpha = prior_alpha,
#'   a_r     = prior_r
#' )
#'
#' ## If all prior values are the same, they can be specified in shorthand
#' ## Equivalent result to the longform priors specified above
#' priors <- list(
#'   a_theta = 0.01,
#'   b_theta = 0.00001,
#'   a_alpha = 1,
#'   a_r     = 0.1
#' )
#'
#' #### Set initial values (optional) ##########################################
#' types  <- unique(sim_SA$cases$Type)
#' q_long <- data.frame(q=rep(15, length(types)), Type=types)
#' init_q <- Q(q_long, q = 'q', type = 'Type')
#' inits <- list(q = init_q) # Pack starting values into a list
#'
#' #### Construct model ########################################################
#' my_model <- HaldDP(y = y, x = x, k = k, priors = priors, inits = inits, a_q = 0.1)
#'
#' #### Set mcmc parameters ####################################################
#' my_model$mcmc_params(n_iter = 2, burn_in = 2, thin = 1)
#'
#' #### Update model ###########################################################
#' my_model$update()
#' ## Add an additional 10 iterations
#' my_model$update(n_iter = 2, append = TRUE)
#'
#' #### Extract posterior ######################################################
#' ## returns the posterior for the r, alpha, q, c,
#' ## lambda_i, xi and xi_prop parameters,
#' ## for all times, locations, sources and types
#' ## the posterior is returned as a list or arrays
#' \dontrun{str(my_model$extract())}
#'
#' ## returns the posterior for the r and alpha parameters,
#' ## for time 1, location B, sources Source3, and Source4,
#' ## types 5, 25, and 50, and iterations 200:300
#' ## the posterior is returned as a list of dataframes
#' \dontrun{
#' str(my_model$extract(params = c("r", "alpha"),
#'                  times = "1", location = "B",
#'                  sources = c("Source3", "Source4"),
#'                  types = c("5", "25", "50"),
#'                  iters = 5:15,
#'                  flatten = TRUE))
#' }
#'
#' #### Calculate medians and credible intervals ###############################
#' \dontrun{my_model$summary(alpha = 0.05, CI_type = "chen-shao")}
#' ## subsetting is done in the same way as extract()
#' \dontrun{my_model$summary(alpha = 0.05, CI_type = "chen-shao",
#'                  params = c("r", "alpha"),
#'                  times = "1", location = "B",
#'                  sources = c("Source3", "Source4"),
#'                  types = c("5", "25", "50"),
#'                  iters = 5:15,
#'                  flatten = TRUE)
#' }
#'
#' #### Plot heatmap and dendrogram of the type effect grouping ################
#' my_model$plot_heatmap()
#'
#' #### Extract data, initial values, prior values, acceptance
#' ## rates for the mcmc algorithm and mcmc parameters
#' my_model$get_data()
#' my_model$get_inits()
#' my_model$get_priors()
#' my_model$get_acceptance()
#' my_model$get_mcmc_params()
#'
#'
#'

HaldDP = function(y, x, k, priors, a_q, inits = NULL)
  HaldDP_$new(y, x, k, priors, a_q, inits)

HaldDP_ <- R6::R6Class(
  "HaldDP",
  private = list(
    ## DATA
    y = NULL,
    X = NULL,
    ## 3D array [sources, type, time] giving the number of positive samples
    k = NULL,

    nTypes = NULL,
    nSources = NULL,
    nTimes = NULL,
    nLocations = NULL,

    namesTypes = NULL,
    namesSources = NULL,
    namesTimes = NULL,
    namesLocations = NULL,

    priors = NULL,
    inits = NULL,
    # Chain starting values

    ## MCMC PARAMETERS
    n_iter = NULL,
    n_iter_old = NULL,
    append = NULL,
    burn_in = NULL,
    thin = NULL,
    total_iters = 0,
    n_r = NULL,
    update_schema = NULL,

    ## MODEL and RESULTS
    modelParams = c('alpha', 'r', 'q'),
    a_q = NULL,
    DPModel_impl = NULL,
    posterior = NULL,
    acceptance = NULL,

    ## MCMC updaters for each parameter
    updaters = NULL,

    # Methods
    reset_posterior = function()
    {
      private$total_iters = 0
      private$posterior = Posterior_HaldDP$new(
        n_iter = private$n_iter,
        namesSources = private$namesSources,
        namesTimes = private$namesTimes,
        namesLocations = private$namesLocations,
        namesTypes = private$namesTypes
      )
    },
    set_update_schema = function(update_schema)
    {
      if (!all(update_schema %in% private$modelParams))
        stop("Unknown parameter in update_schema")
      if (!setequal(private$update_schema, update_schema) &
          isTRUE(private$append))
        stop("Append must be set to false if changing update_schema, model parameters, or data.")
      private$update_schema = update_schema
      private$assign_updaters()
    },
    assign_updaters = function()
    {
      steps = list()
      for (time in 1:private$nTimes) {
        for (location in 1:private$nLocations) {
          if ('alpha' %in% private$update_schema) {
            steps[[length(steps) + 1]] =
              AdaptiveLogDirMRW$new(
                private$DPModel_impl$alphaNodes[[time]][[location]],
                tune = rep(0.4, private$nSources),
                name = paste0(
                  "alpha: time ",
                  private$namesTimes[time],
                  ", location ",
                  private$namesLocations[location]
                )
              ) # TODO: User specificed starting variance?
          }
        }

        if ('r' %in% private$update_schema) {
          for (source in 1:private$nSources) {

            # method of moments: beta to choose tuning values
            alpha <-
              private$priors$a_r[, source, time] + private$X[, source, time]
            alpha_0 <- sum(alpha)
            var_alphas <-
              (alpha * (alpha_0 - alpha)) / ((alpha_0 ^ 2) * (alpha_0 + 1))
            tune_val <- 100 * sqrt(var_alphas)
            steps[[length(steps) + 1]] =
              AdaptiveLogDirMRW$new(
                private$DPModel_impl$rNodes[[time]][[source]],
                toupdate = function()
                  sample(private$nTypes, private$n_r),
                tune = tune_val,
                batchsize = 10,
                name = paste0(
                  "r: time ",
                  private$namesTimes[time],
                  ", source ",
                  private$namesSources[source]
                )
              )
          }
        }
      }
      if ('q' %in% private$update_schema) {
        steps[[length(steps) + 1]] <-
          PoisGammaDPUpdate$new(private$DPModel_impl$qNodes)
      }
      private$updaters = steps
    },
    checkin_data = function(y,x,k)
    {
      # Check dimension mappings
      if(!identical(dimnames(y$x)$type, dimnames(x$x)$type))
        stop("Types in x and y do not match")
      if(!identical(dimnames(y$x)$time, dimnames(x$x)$time))
        stop("Times in x and y do not match")
      if(!identical(dimnames(x$x)$time, dimnames(k$x)$time))
        stop("Times in x and k do not match")
      if(!identical(dimnames(x$x)$source, dimnames(k$x)$source))
        stop("Sources in x and k do not match")

      private$y = y$x
      private$X = x$x
      private$k = k$x

      private$nTypes     = dim(private$y)[1]
      private$nTimes     = dim(private$y)[2]
      private$nLocations = dim(private$y)[3]
      private$nSources   = dim(private$X)[2]

      private$namesTypes     = dimnames(private$y)$type
      private$namesTimes     = dimnames(private$y)$time
      private$namesLocations = dimnames(private$y)$location
      private$namesSources   = dimnames(private$X)$source
    },
    checkin_a_q = function(a_q)
    {
      if (length(a_q) > 1)
        warning("length(a_q) > 1 so only first element will be used")
      if (!isFinitePositive(a_q[1]) |
          a_q[1] <= 0)
        stop("a_q should be a single positive number")
      private$a_q <- a_q[1]
    },
    checkin_priors = function(priors)
    {
      if (!(is.list(priors) &
          all(c("a_alpha", "a_r", "a_theta", "b_theta") %in% names(priors))))
        stop("The priors must be a list with names a_alpha, a_r, a_theta, and b_theta.")

      ## theta priors
      if (length(priors$a_theta) != 1 |
          !isFinitePositive(priors$a_theta))
        stop("priors$a_theta must be a single positive number.")
      if (length(priors$b_theta) != 1 |
          !isFinitePositive(priors$b_theta))
        stop("priors$b_theta must be a single positive number.")

      ## r priors
      if(is.atomic(priors$a_r) & length(priors$a_r)==1) {
        if(!isFinitePositive(priors$a_r))
          stop("priors$a_r must be positive for a single number specification.")
        priors$a_r <- array(
          priors$a_r,
          dim = c(private$nTypes, private$nSources, private$nTimes),
          dimnames = list(
            type = private$namesTypes,
            source = private$namesSources,
            time = private$namesTimes
          )
        )
      }
      else if(identical(class(priors$a_r), c('X','Data','R6')))
      {
        if(!identical(dimnames(priors$a_r$x), dimnames(private$X)))
          stop("dimnames(priors$a_r$x) must match dimnames(X)")
        priors$a_r = priors$a_r$x
      }
      else
        stop("Prior for r must be either a single number or of type sourceR::X")

      ## alpha priors
      if(is.atomic(priors$a_alpha) & length(priors$a_alpha)==1) {
        if(!isFinitePositive(priors$a_alpha))
          stop("priors$a_alpha <= 0.  Must be positive.")
        priors$a_alpha = array(
          priors$a_alpha,
          dim = c(private$nSources, private$nTimes, private$nLocations),
          dimnames = list(
            source = private$namesSources,
            time = private$namesTimes,
            location = private$namesLocations
          )
        )
      }
      else if(identical(class(priors$a_alpha), c('Alpha','Data','R6')))
      {
        dn = list(source=private$namesSources, time=private$namesTimes, location=private$namesLocations)
        if(!identical(dimnames(priors$a_alpha$x), dn))
          stop("priors$a_alpha must have source/time/location dimensions appropriate to data.")
         priors$a_alpha = priors$a_alpha$x
      }
      else
        stop("Prior parameters for alpha must be specified either as a constant or of type sourceR::Alpha")
      # Stow as a member
      private$priors = priors
    },
    checkin_inits = function(inits)
    {
      private$inits = list()
      ## r values
      if (!("r" %in% names(inits))) {
        ## default is the source data matrix MLE (plus stability factor)
        private$inits$r = apply(private$X + 1e-6, c('source', 'time'), function(x)
          x / sum(x)) %>% aperm(c('type', 'source', 'time'))
      } else {
        if (!identical(class(inits$r), c('X','Data','R6')))
          stop("inits$r must be an object of type sourceR::X")
        if (!identical(dimnames(inits$r$x), dimnames(private$X)))
          stop("inits$r must be of the same rank as source data")
        if (any(inits$r$x <=0))
          stop("One or more elements of initial r values <= 0.")
        srcSums = apply(inits$r$x, c('source','time'), sum)
        if(any(srcSums != 1))
          warning("1 or more r[,source,time] vectors did not sum to 1.  Renormalising.")
        private$inits$r = apply(inits$r$x, c('source','time'), function(x) x/sum(x))
      }
      ## alpha inits
      if (!("alpha" %in% names(inits))) {
          private$inits$alpha = apply(private$priors$a_alpha, c('time', 'location'), function(x) {
          z = gtools::rdirichlet(1, x)
          names(z) = names(x)
          z
        })
      } else {
        if(!identical(class(inits$alpha), c('Alpha','Data','R6')))
          stop('inits$alpha must be of class sourceR::Alpha')
        if(!( identical(dimnames(inits$alpha$x)$source, private$namesSources) &
              identical(dimnames(inits$alpha$x)$time, private$namesTimes) &
              identical(dimnames(inits$alpha$x)$location, private$namesLocations)))
          stop('inits$alpha must have same number of sources and times as X, and same number of locations as y')
        srcSums = apply(inits$alpha$x, c('time','location'), sum)
        if(any(srcSums != 1))
          warning("1 or more init$alpha vectors did not sum to 1. Renormalising.")
        private$inits$alpha = apply(inits$alpha$x, c('time','location'), function(x) x/sum(x))
      }

      ## q inits
      if (!("q" %in% names(inits))) {
        ## default to all in 1 group with mean solved using the initial parameters for r and alpha
        ## sum_{i=1}^n lambda_i = q * sum_{t=1}^T sum_{l=1}^L (sum_{i=1}^n sum_{j=1}^m p_{ijt} a_{jtl})
        ## sum_{i=1}^n y_i = q * sum_{t=1}^T sum_{l=1}^L (sum_{i=1}^n sum_{j=1}^m p_{ijt} a_{jtl})
        ## q = sum_{i=1}^n y_i / sum_{t=1}^T sum_{l=1}^L (sum_{i=1}^n sum_{j=1}^m p_{ijt} a_{jtl})
        tt = private$inits[c('alpha','r')]
        tt$k = private$k
        tt = lapply(tt, function(x) as.tensor.default(x, dims=dimnames(x)))
        lambdaBar = mul.tensor(
          tt$r,
          i = 'source',
          tt$alpha * tt$k,
          j = 'source',
          by = c('type','time','location')
        ) %>% sum
        q_val <- sum(private$y) / lambdaBar
        private$inits$theta <- q_val
        private$inits$s <- rep(1, private$nTypes)
        names(private$inits$s) = private$namesTypes
      } else {
        if(!identical(class(inits$q), c('Q','Data','R6')))
          stop("q initialiser must be of class sourceR::Q")
        if(!identical(names(inits$q$s), private$namesTypes))
          stop("q initialiser must have type names identical to y")
        private$inits$theta = inits$q$theta
        private$inits$s = inits$q$s
      }
    },
    set_niter = function(n_iter)
    {
      assert_that(is.numeric(n_iter), n_iter > 0, is.finite(n_iter))
      private$n_iter = n_iter
    },
    set_append = function(append)
    {
      assert_that(is.logical(append), is.atomic(append))
      private$append <- append
    },
    set_burn_in = function(burn_in)
    {
      if (!isFiniteInteger(burn_in) | burn_in < 0) {
        stop("burn_in is not a positive integer.")
      } else {
        private$burn_in <- burn_in
      }
    },
    set_thin = function(thin)
    {
      assert_that(isFinitePositiveInteger(thin))
      private$thin <- thin
    },
    set_n_r = function(n_r)
    {
      assert_that(isFinitePositiveInteger(n_r),
                  n_r <= private$nTypes)
      private$n_r <- n_r
    },
    calc_CI = function(x, alpha, CI_type)
    {
      # x is the MCMC output
      switch(
        CI_type,
        "chen-shao" = ci_chenShao(x, alpha),
        "percentiles" = ci_percentiles(x, alpha),
        "SPIn" = ci_SPIn(x, alpha),
        stop(
          "The type of interval specified must be one of: chen-shao, percentiles, or SPIn."
        )
      )
    },
    calc_CI_alpha = function(object, alpha, CI_type)
    {
      res = as.tensor(apply(object, c(1,2,3),
                            private$calc_CI,
                            alpha=alpha,
                            CI_type=CI_type))
      names(res)[1] = 'ci'
      dimnames(res)$ci = c('lower','median','upper')
      res
    },
    calc_CI_q = function(object, alpha, CI_type)
    {
      res = as.tensor(apply(object, 'type', private$calc_CI, alpha=alpha, CI_type=CI_type))
      names(res)[1] = 'ci'
      dimnames(res)$ci = c('lower','median','upper')
      res
    },
    calc_CI_r = function(object, alpha, CI_type)
    {
      res = as.tensor(apply(object, c('type','source','time'),
                            private$calc_CI,
                            alpha=alpha,
                            CI_type=CI_type))
      names(res)[1] = 'ci'
      dimnames(res)$ci = c('lower','median','upper')
      res
    },
    calc_CI_lambda_i = function(object, alpha, CI_type)
    {
      res = to.tensor(apply(object, c('type','time','location'),
                            private$calc_CI,
                            alpha=alpha,
                            CI_type=CI_type))
      names(res)[1] = 'ci'
      dimnames(res)$ci = c('lower','median','upper')
      res
    },
    calc_CI_xi = function(object, alpha, CI_type)
    {
      res = as.tensor(apply(object, c(1,2,3),
                            private$calc_CI,
                            alpha=0.05,
                            CI_type=CI_type))
      names(res)[1] = 'ci'
      dimnames(res)$ci = c('lower','median','upper')
      res
    },
    check_extract_summary_params = function(params,
                                            times,
                                            locations,
                                            sources,
                                            types,
                                            iters,
                                            flatten)
    {
      if (!mode(params) %in% c("character") |
          length(params) > 7 |
          length(unique(params)) != length(params) |
          !all(
            unique(params) %in% c(
              "alpha",
              "q",
              "s",
              "r",
              "lambda_i",
              "xi",
              "xi_prop"
            )
          )) {
        stop(
          paste(
            "params must be a vector where elements are a subset of ",
            paste(
              c(
                "alpha",
                "q",
                "s",
                "r",
                "lambda_i",
                "xi",
                "xi_prop"
              ),
              collapse = ", "
            ),
            " with no repeated numbers."
          )
        )
      }

      if (is.null(iters)) {
        ## return all iters
        iters <- 1:private$posterior$iters
      } else {
        if (!mode(iters) %in% c("numeric") |
            !all(isFinitePositiveInteger(iters)) |
            !all(iters <= private$posterior$iters) |
            !length(unique(iters)) == length(iters)) {
          stop(
            "iters must be a vector where all elements are integers between 0 and n_iter, with no repeated numbers."
          )
        }
      }

      if (is.null(times)) {
        ## return all times
        times <- private$namesTimes
      } else {
        if (!mode(times) %in% c("numeric", "character") |
            length(times) > private$nTimes |
            length(unique(times)) != length(times) |
            !all(unique(times) %in% private$namesTimes)) {
          stop(
            paste(
              "times must be a vector where elements can only be ",
              paste(private$namesTimes, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (is.null(locations)) {
        ## return all locations
        locations <- private$namesLocations
      } else {
        if (!mode(locations) %in% c("numeric", "character") |
            length(locations) > private$nLocations |
            length(unique(locations)) != length(locations) |
            !all(unique(locations) %in% private$namesLocations)) {
          stop(
            paste(
              "locations must be a vector where elements can only be ",
              paste(private$namesLocations, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (is.null(sources)) {
        ## return all sources
        sources <- private$namesSources
      } else {
        if (!mode(sources) %in% c("numeric", "character") |
            length(sources) > private$nSources |
            length(unique(sources)) != length(sources) |
            !all(unique(sources) %in% private$namesSources)) {
          stop(
            paste(
              "sources must be a vector where elements can only be ",
              paste(private$namesSources, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (is.null(types)) {
        ## return all sources
        types <- private$namesTypes
      } else {
        if (!mode(types) %in% c("numeric", "character") |
            length(types) > private$nTypes |
            length(unique(types)) != length(types) |
            !all(unique(types) %in% private$namesTypes)) {
          stop(
            paste(
              "types must be a vector where elements can only be ",
              paste(private$namesTypes, collapse = ", "),
              " with no repeated values."
            )
          )
        }
      }

      if (length(flatten) != 1 |
          !isFiniteLogical(flatten))
        stop("flatten must be a single logical value.")

      ## calculate the lambda's
      if ("lambda_i" %in% params) {
        private$posterior$calc_lambda_i(
          private$posterior$iters,
          private$nTimes,
          private$nLocations,
          private$nTypes,
          private$namesTimes,
          private$namesLocations,
          private$namesTypes,
          1:private$posterior$iters,
          private$k
        )
      }
      if ("xi" %in% params) {
        private$posterior$calc_xi(
          private$posterior$iters,
          private$nSources,
          private$nTimes,
          private$nLocations,
          private$namesSources,
          private$namesTimes,
          private$namesLocations,
          1:private$posterior$iters,
          private$k
        )
      }

      if ("xi_prop" %in% params) {
        private$posterior$calc_xi_prop(
          private$posterior$iters,
          private$nSources,
          private$nTimes,
          private$nLocations,
          private$namesSources,
          private$namesTimes,
          private$namesLocations,
          1:private$posterior$iters,
          private$k
        )
      }

      return(
        list(
          params = params,
          times = times,
          locations = locations,
          sources = sources,
          types = types,
          iters = iters,
          flatten = flatten
        )
      )
    }
  ),
  public = list(
    initialize = function(y, x, k, priors, a_q, inits)
    {
      if(!identical(class(y), c("Y","Data","R6")) )
        stop("y must be of type sourceR::Y")
      if(!identical(class(x), c("X","Data","R6")) )
        stop("x must be of type sourceR::X")
      if(!identical(class(x), c("Prev","Data","R6")))
      private$checkin_data(y,x,k)

      ## read in priors
      private$checkin_priors(priors)
      private$checkin_a_q(a_q)

      ## read in initial values (must go after priors as uses priors to generate initial values)
      private$checkin_inits(inits)

      private$DPModel_impl =
        DPModel_impl$new(
          y = private$y,
          X = private$X,
          R = private$inits$r,
          Time = private$namesTimes,
          Location = private$namesLocations,
          Sources = private$namesSources,
          Type = private$namesTypes,
          prev = private$k,
          a_q = private$a_q,
          a_theta = private$priors$a_theta,
          b_theta = private$priors$b_theta,
          a_r = private$priors$a_r,
          a_alpha = private$priors$a_alpha,
          s = private$inits$s,
          theta = private$inits$theta,
          alpha = private$inits$alpha
        )
      self$mcmc_params()
    },
    mcmc_params = function(n_iter = 1000,
                          burn_in = 0,
                          thin = 1,
                          n_r = ceiling(private$nTypes * 0.2),
                          update_schema = c('q', 'alpha', 'r'))
    {
      private$set_append(FALSE)
      private$set_niter(n_iter)
      private$set_burn_in(burn_in)
      private$set_thin(thin)
      private$set_n_r(n_r)
      private$set_update_schema(update_schema) # Builds the schema here
      private$reset_posterior() # Resets posterior counters
    },
    # Run the MCMC.  If n_iter is specified,
    # runs the algorithm until n_iter iterations
    # have been added to the posterior, taking into
    # account thinning and burn-in.  Adds to an existing
    # posterior if append is TRUE.
    update = function(n_iter, append = TRUE)
    {
      iters = 0
      if (!missing(n_iter)) {
        assert_that(isFinitePositiveInteger(n_iter))
        iters = n_iter
      }
      else
        iters = private$n_iter

      ## make posterior if append = F, otherwise extend its size
      assert_that(is.logical(append))
      if (!append)
        private$reset_posterior()

      ## run mcmc
      ## user feedback and main MCMC loop
      pb <- txtProgressBar(max = iters, style = 3)
      i = 0
      while (i < iters) {
        ## save results if iteration i is > burn in and a multiple of thin
        if (private$total_iters > private$burn_in &
            private$total_iters %% private$thin == 0)
        {
          private$posterior$sample(private$DPModel_impl)
          i = i + 1
        }
        # Do the updates
        lapply(private$updaters, function(x)
          x$update())

        private$total_iters = private$total_iters + 1
        setTxtProgressBar(pb, i)
      }
    },

    ## Functions to access the data and results
    get_data = function()
    {
        list(
          y = private$y,
          X = private$X,
          k = private$k
        )
    },
    get_priors = function()
    {
     list(priors = private$priors,
                  a_q = private$a_q)
    },
    get_inits = function()
    {
      private$inits
    },
    get_mcmc_params = function()
    {
        list(
          n_iter = private$n_iter,
          append = private$append,
          burn_in = private$burn_in,
          thin = private$thin,
          update_schema = private$update_schema
        )
    },
    get_acceptance = function()
    {

      ## mcmc is finished, save acceptance rate summary for r and alpha
      acceptance =
        Acceptance$new(
          nSources = private$nSources,
          nTimes = private$nTimes,
          nLocations = private$nLocations,
          nTypes = private$nTypes,
          namesSources = private$namesSources,
          namesTimes = private$namesTimes,
          namesLocations = private$namesLocations,
          namesTypes = private$namesTypes,
          updateSchema = private$update_schema
        )

      sapply(1:length(private$updaters), function(i) {
        tryCatch({
          acceptances <- private$updaters[[i]]$acceptance()
          names <- private$updaters[[i]]$name

          tmp <- strsplit(names, split = " ")[[1]]
          if (tmp[1] == "r:") {
            t_name <- substr(tmp[3], 1, nchar(tmp[3]) - 1)
            s_name <- tmp[5]
            acceptance$r[,
                         which(s_name == dimnames(acceptance$r)$source),
                         which(t_name == dimnames(acceptance$r)$time)] <-
              acceptances
          } else if (tmp[1] == "alpha:") {
            t_name <- substr(tmp[3], 1, nchar(tmp[3]) - 1)
            l_name <- tmp[5]
            acceptance$alpha[,
                             which(t_name == dimnames(acceptance$alpha)$time),
                             which(l_name == dimnames(acceptance$alpha)$location)] <-
              acceptances
          }

        },
        error = function(e) {
          NULL
        })
      })
      list(alpha=acceptance$alpha, r=acceptance$r)
    },
    extract = function(params = c("alpha", "q", "s", "r", "lambda_i", "xi", "xi_prop"),
                       times = NULL,
                       locations = NULL,
                       sources = NULL,
                       types = NULL,
                       iters = NULL,
                       flatten = FALSE,
                       drop = TRUE)
    {
      params_checked =
        private$check_extract_summary_params(params, times, locations,
                                             sources, types, iters, flatten)

      assert_that(is.atomic(drop), is.logical(drop))
      if(private$posterior$iters == 0) warning("Empty posterior")
      params = params_checked$params
      times = params_checked$times
      locations = params_checked$locations
      sources = params_checked$sources
      types = params_checked$types
      iters = params_checked$iters
      flatten = params_checked$flatten

      if(isTRUE(flatten)) drop = FALSE
      res = lapply(setNames(params, params), function(x) {
        switch(
          x,
          # TODO: This is a PITA hack for buggy tensorA::[.tensor
          # method: replaced with sliceTensor(x,...) for now.
          # See utils.R for details. CPJ 2017-04-06
          "alpha" = sliceTensor(
            private$posterior$alpha,
            sources,
            times,
            locations,
            iters,
            drop = drop
          ),
          "q"     = sliceTensor(private$posterior$q, types, iters, drop = drop),
          "s"     = sliceTensor(private$posterior$s, types, iters, drop = drop),
          "r"     = sliceTensor(private$posterior$r, types, sources, times, iters, drop = drop),
          "lambda_i" = sliceTensor(
            private$posterior$lambda_i,
            types,
            times,
            locations,
            iters,
            drop = drop
          ),
          "xi" = sliceTensor(
            private$posterior$xi,
            sources,
            times,
            locations,
            iters,
            drop = drop
          ),
          "xi_prop" = sliceTensor(
            private$posterior$xi_prop,
            sources,
            times,
            locations,
            iters,
            drop = drop
          ),
          stop("Unrecognised model component")
        )
      })

      if (isTRUE(flatten)) {
        res = lapply(res, reshape2::melt)
      }
      res
    },
    plot_heatmap = function(iters,
                            cols = c("blue", "white"),
                            hclust_method = "complete")
    {
      if (! ("q" %in% private$update_schema)) {
        stop(
          "q is not in private$update_schema, therefore q was not updated nd all values of the marginal posterior for q and s are the same (and a heatmap can't be plotted)."
        )
      }
      if (missing(iters)) {
        ## return all iters
        iters <- 1:private$posterior$iters
      } else {
        if (!mode(iters) %in% c("numeric") |
            !all(isFinitePositive(iters)) |
            !all(isFiniteInteger(iters)) |
            !all(iters <= private$posterior$iters) |
            !length(unique(iters)) == length(iters)) {
          stop(
            "iters must be a vector where all elements are integers between 0 and n_iter, with no repeated numbers."
          )
        }
      }

      # Draw the heatmap
      groups <- as.data.frame(private$posterior$s[, iters], stringsAsFactors = T)
      clusterHeatMap(groups, cols, rownames(private$posterior$q), hclust_method)
    },
    summary = function(alpha = 0.05,
                       params = c("alpha", "q", "r", "lambda_i", "xi", "xi_prop"),
                       times = NULL,
                       locations = NULL,
                       sources = NULL,
                       types = NULL,
                       iters = NULL,
                       flatten = FALSE,
                       CI_type = "chen-shao")

    {
      if(private$posterior$iters == 0) stop("Cannot summarise empty posterior.
                                            Please call update() first.")
      object =
        self$extract(params, times, locations, sources, types, iters, FALSE, drop = F)

      s = lapply(setNames(params, params), function(x) {
        switch(
          x,
          "alpha" = private$calc_CI_alpha(object$alpha, alpha, CI_type),
          "q"     = private$calc_CI_q(object$q, alpha, CI_type),
          "r"     = private$calc_CI_r(object$r, alpha, CI_type),
          "lambda_i" = private$calc_CI_lambda_i(object$lambda_i, alpha, CI_type),
          "xi" = private$calc_CI_xi(object$xi, alpha, CI_type),
          "xi_prop" = private$calc_CI_xi(object$xi_prop, alpha, CI_type)
        )
      })

      if (flatten == TRUE) {
        s = lapply(s, function(x) {
          vars = names(x)[names(x)!='ci']
          frm = paste(vars, collapse='+') %>% paste0('~ci')
          class(x) = 'array'
          reshape2::melt(x) %>% reshape2::dcast(eval(parse(text=frm)), value.var='value')
        })
      }
      s
    }
  )
)
