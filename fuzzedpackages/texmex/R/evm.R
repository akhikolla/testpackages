#' Extreme value modelling
#'
#' Likelihood based modelling and inference for extreme value models, possibly
#' with explanatory variables.
#'
#' The main modelling function is \code{evm} (extreme value model) and the
#' distribution to be used is specified by passing an object of class
#' \code{texmexFamily} to the \code{family} argument.
#'
#' The default \code{texmexFamily} object used by \code{evm} is \code{gpd}.
#' Currently, the other \code{texmexFamily} objects available are \code{gev}
#' which results in fitting a generalized extreme value (GEV) distribution to
#' the data, \code{gpdIntCensored} which can be used to fit the GPD to data which has
#' been rounded to a given numebr of decimal places by recognisiing the data as
#' interval censored, and \code{egp3} which fits the extended generalized Pareto
#' distribution version 3 of Papastathopoulos and Tawn (2013).
#'
#' See Coles (2001) for an introduction to extreme value modelling and the GPD
#' and GEV models.
#'
#' For the GPD model, we use the following parameterisation of evm:
#'
#' \deqn{P(Y \le y) = 1 - (1 + \xi y / \sigma)^{-1/\xi}}
#' for \eqn{y \ge 0} and \eqn{1 + \xi y / \sigma \ge 0.}
#'
#' For the GEV model, we use:
#'
#' \deqn{P(Y \le y) = exp (-(1 + \xi (y - \mu) / \sigma) ^{-1/\xi})}
#'
#'
#' In each case, the scale parameter is sigma (\eqn{\sigma}) and the shape
#' parameter is xi (\eqn{\xi}). The GEV distribution also has location
#' parameter mu (\eqn{\mu}). See Papastathopoulos and Tawn (2013) for
#' specification of the EGP3 model.
#'
#' Working with the log of the scale parameter improves the stability of
#' computations, makes a quadratic penalty more appropriate and enables the
#' inclusion of covariates in the model for the scale parameter, which must
#' remain positive.  We therefore work with \eqn{\phi}=log(\eqn{\sigma}).  All
#' specification of priors or penalty functions refer to \eqn{\phi} rather than
#' \eqn{\sigma}.  A quadratic penalty can be thought of as a Gaussian prior
#' distribution, whence the terminology of the function.
#'
#' Parameters of the evm are estimated by using maximum (penalized) likelihood
#' (\code{method = "optimize"}), or by simulating from the posterior
#' distribution of the model parameters using a Metropolis algorithm
#' (\code{method = "simulate"}).  In the latter case, \code{start} is used as a
#' starting value for the Metropolis algorithm; in its absence, the maximum
#' penalized likelhood point estimates are computed and used.
#'
#' A boostrap approach is also available (\code{method = "bootstrap"}).  This
#' runs a parametric bootstrap, simulating from the model fit by optimization.
#'
#' When \code{method = "simulate"} the \code{print} and \code{summary}
#' functions give posterior means and standard deviations. Posterior means are
#' also returned by the \code{coef} method. Depending on what you want to do
#' and what the posterior distributions look like (use \code{plot} method) you
#' might want to work with quantiles of the posterior distributions instead of
#' relying on standard errors.
#'
#' When \code{method = "bootstrap"}, summaries of the bootstrap distribution
#' and the bootstrap estimate of bias are displayed.
#'
#' @param y Either a numeric vector or the name of a variable in \code{data}.
#' @param data A data frame containing \code{y} and any covariates.
#' @param family An object of class 'texmexFamily'. Defaults to
#' \code{family=gpd} and a generalized Pareto distribution (GPD) is fit to the data.
#' Alternatively the family could be \code{gev}, \code{weibull} or
#' \code{gumbel}, resulting in a generalized extreme value distribution, Weibull
#' or Gumbell distribution being fit. Family \code{cgpd} fits the generalized
#' Pareto distribution but with the shape parameter constrained to be
#' > 0.5 by using the link function suggested by Yee and Stephenson (2007),
#' \eqn{\eta} = log(\eqn{\xi} + 0.5). Family \code{egp3} fits the extended
#' GP family 3 of Papastathopoulos and Tawn (2013). No other families are currently
#' available in texmex, but users may write their own.
#' @param th For threshold excess models (such as when \code{family=gpd}), the
#' threshold for \code{y}, exceedances above which will be used to fit the
#' upper tail model. Note that if you have already thresholded your data and
#' want to model all of \code{y}, you still need to specify \code{th}.
#' @param qu An alternative to \code{th}, a probability defined such that
#' \code{quantile(y,qu)} equals \code{th}.
#' @param ...  In \code{evm}, formulae for the parameters in the family, e.g.
#' \code{phi = ~ x}. If none are specified, they all default to \code{~1}.
#' @param penalty How to penalize the likelhood. Currently, either "none"",
#' "gaussian"" or "lasso"" are the only allowed values. If \code{penalty} is
#' "gaussian" or "lasso" then the parameters for the penalization are specified
#' through the \code{priorParameters} argument. See below. Defaults to
#' \code{penalty=NULL} and applies maximum likelihood estimation.
#' @param prior If \code{method = "optimize"}, just an alternative way of
#' specifying the penalty, and only one or neither of \code{penalty} and
#' \code{prior} should be given. If \code{method = "simulate"}, prior must be
#' "gaussian" because no other prior distributions have been implemented.
#' @param method Should be either "optimize" (the default), "simulate" or
#' "bootstrap".  The first letter or various abbreviations will do. If
#' 'optimize' is used, the (penalized) likelihood is directly optimized using
#' \code{optim} and point estimates (either ML or MAP estimates) are returned
#' with other information. If "simulate", a Metropolis algorithm is used to
#' simulate from the joint posterior distribution of the parameters. If
#' "bootstrap", a parametric boostrap is performed.
#' @param cov How to compute the covariance matrix of the parameters. Defaults
#' to \code{cov = "observed"} in which case the observed information matrix is
#' used, if the \code{info} element of the \code{texmexFamily} object is
#' present. Note that currently, this is not implemented for \code{gev}.
#' Alternatives are \code{cov = "numeric"} in which case a numerical
#' approximation of the Hessian is used (see the help for \code{optim}), or
#' \code{cov = "sandwich"} if the \code{sandwich} element of the
#' \code{texmexFamily} object is implemented.  The \code{cov = "sandwich"}
#' method implements the Huber sandwich correction to the covariance matrix for
#' data which are not independent and in which case the likelihood function no
#' longer has the interpretation of a joint likelihood, but instead should be
#' interpreted as a pseudo-likelihod.
#'
#' In some cases, particularly with small samples, the numerical approximation
#' can be quite different from the closed form (\code{cov="observed"}) result,
#' and the value derived from the observed information should be preferred.
#' However, in either case, since the underlying log-likelihood may be far from
#' quadratic for small samples, the resulting estimates of standard errors are
#' liable to approximate poorly the true standard errors. Also see the comments
#' in the Details section, below.
#' @param start Starting values for the parameters, to be passed to
#' \code{optim}.  If not provided, the function will use the \code{start}
#' element of the \code{texmexFamily} object if it exists.
#' @param priorParameters A list with two components. The first should be a
#' vector of means, the second should be a covariance matrix if the
#' penalty/prior is "gaussian" or "quadratic" and a diagonal precision matrix
#' if the penalty/prior is "lasso", "L1" or "Laplace".  If \code{method =
#' "simulate"} then these represent the parameters in the Gaussian prior
#' distribution.  If \code{method = 'optimize'} then these represent the
#' parameters in the penalty function.  If not supplied: all default prior
#' means are zero; all default prior variances are \eqn{10^4}; all covariances
#' are zero.
#' @param maxit The number of iterations allowed in \code{optim}.
#' @param trace Whether or not to print progress to screen. If \code{method =
#' "optimize"}, the argument is passed into \code{optim} -- see the help for
#' that function. If \code{method = "simulate"}, the argument determines at how
#' many steps of the Markov chain the function should tell the user, and in
#' this case it defaults to \code{trace = 10000}.
#' @param iter Number of simulations to generate under \code{method =
#' "simulate"}.  Defaults to 40500.
#' @param burn The number of initial steps to be discarded. Defaults to 500.
#' @param thin The degree of thinning of the resulting Markov chains. Defaults
#' to 4 (one in every 4 steps is retained).
#' @param chains The number of Markov chains to run. Defaults to 1. If you run
#'   more than 1, the function tries to figure out how to do it in parallel
#'   using as many cores as there are chains.
#' @param proposal.dist The proposal distribution to use, either multivariate
#' gaussian or a multivariate Cauchy.
#' @param jump.cov Covariance matrix for proposal distribution of Metropolis
#' algorithm.  This is scaled by \code{jump.const}.
#' @param jump.const Control parameter for the Metropolis algorithm.
#' @param verbose Whether or not to print progress to screen. Defaults to
#' \code{verbose=TRUE}.
#' @param R The number of parametric bootstrap samples to run when \code{method
#' = "bootstrap"} is requested. Defaults to 1000.
#' @param cores The number of cores to use when bootstrapping. Defaults to
#' \code{cores=NULL} and the function guesses how many cores are available and
#' uses them all.
#' @param export Character vector of names of objects to export if parallel
#'   processing is being used and you are using objects from outside of
#'   texmex. It it passed to \code{parallel::clusterExport} and used by
#'   \code{texmex::evmBoot}.
#' @return If \code{method = "optimize"}, an object of class \code{evmOpt}:
#'
#' \item{call}{The call to \code{evmSim} that produced the object.}
#' \item{data}{The original data (above and below the threshold for fitting if
#' a distribution for threshold excesses has been used). In detail, \code{data}
#' is a list with elements \code{y} and \code{D}. \code{y} is the response
#' variable and \code{D} is a list containing the design matrices implied by
#' any formlae used in the call to \code{evm}.} \item{convergence}{Output from
#' \code{optim} relating to whether or not the optimizer converged.}
#' \item{message}{A message telling the user whether or not convergence was
#' achieved.} \item{threshold}{The threshold of the data above which the evmSim
#' model was fit.} \item{penalty}{The type of penalty function used, if any.}
#' \item{coefficients}{The parameter estimates as computed under maximum
#' likelihood or maximum penalized likelihood.} \item{rate}{The proportion of
#' observations above the threshold. If the model is not a threshold exceedance
#' model (e.g. the GEV model), the rate will be 1.}
#'
#' \item{priorParameters}{See above.}
#'
#' \item{residuals}{Residuals computed using the residual function in the
#' \code{texmexFamily} object, if any. These are used primarly for producing
#' QQ and PP plots via \code{\link{plot.evmOpt}} or \code{\link{ggplot.evmOpt}}.
#' The residuals are transformed values of the raw data, accounting for the
#' parameter estimates: see the \code{residuals} component of the
#' \code{\link{texmexFamily}} object for the calculations. For the generalized
#' Pareto family, they are (if the model fits well) standard exponential variates;
#' for the GEV family, standard Gumbel variates.}
#'
#' \item{ploglik}{The value of the
#' optimized penalized log-likelihood.} \item{loglik}{The value of the
#' optimized (unpenalized) log-likelihood. If \code{penalty='none'} is used,
#' this will be identical to \code{ploglik}, above.} \item{cov}{The estimated
#' covariance of the parameters in the model.} \item{se}{The estimated standard
#' errors of the parameters in the model.} \item{xlevels}{A named list
#' containing a named list for each design matrix (main parameter) in the
#' model. Each list contians an element named after each factor in the linear
#' predictor for the respective design matrix. These are used by the
#' \code{predict} method to ensure all factor levels are known, even if they
#' don't appear in \code{newdata}.}
#'
#' If \code{method = "simulate"}, an object of class \code{evmSim}:
#'
#' \item{call}{The call to \code{evmSim} that produced the object.}
#' \item{threshold}{The threshold above which the model was fit.}
#' \item{map}{The point estimates found by maximum penalized likelihood and
#' which were used as the starting point for the Markov chain.  This is of
#' class \code{evmOpt} and methods for this class (such as resid and plot) may
#' be useful.} \item{burn}{The number of steps of the Markov chain that are to
#' be treated as the burn-in and not used in inferences.} \item{thin}{The
#' degree of thinning used.} \item{chains}{The entire Markov chain generated by
#' the Metropolis algorithm.} \item{y}{The response data above the threshold
#' for fitting.} \item{seed}{The seed used by the random number generator.}
#' \item{param}{The remainder of the chain after deleting the burn-in and
#' applying any thinning.}
#'
#' If \code{method = "bootstrap"}, an object of class \code{evmBoot}:
#'
#' \item{call}{The call to \code{evmBoot} that produced the object.}
#' \item{replicates}{The parameter estimates from the bootstrap fits.}
#' \item{map}{The fit by by maximum penalized likelihood to the orginal data.
#' This is of class \code{evmOpt} and methods for this class (such as resid and
#' plot) may be useful.}
#'
#' There are summary, plot, print, residuals and coefficients methods available for these
#' classes.
#' @note For both GPD and GEV models, when there are estimated values of
#' \eqn{\xi \le -0.5}, the regularity conditions of the likelihood break down
#' and inference based on approximate standard errors cannot be performed. In
#' this case, the most fruitful approach to inference appears to be by the
#' bootstrap. It might be possible to simulate from the posterior, but finding
#' a good proposal distribution might be difficult and you should take care to
#' get an acceptance rate that is reasonably high (around 40\% when there are
#' no covariates, lower otherwise). To constrain the parameter space of the GP
#' shape parameter, use \code{family = cgpd} in the call to \code{evm} and
#' the transformation \eqn{\eta} = log(\eqn{\xi} + 0.5) is used, as suggested
#' by Yee and Stephenson (2007).
#' @author Janet E. Heffernan, Harry Southworth. Some of the internal code is
#' based on the \code{gpd.fit} function in the \code{ismev} package and is due
#' to Stuart Coles.
#' @seealso \code{\link{plot.evmOpt}} \code{\link{ggplot.evmOpt}} \code{\link{rl.evmOpt}}, \code{\link{predict.evmOpt}},
#' \code{\link{evm.declustered}}.
#' @references S. Coles. An Introduction to Statistical Modelling of Extreme
#' Values. Springer, 2001.
#'
#' I. Papastathopoulos and J. A. Tawn, Extended generalised Pareto models for
#' tail estimation, Journal of Statistical Planning and Inference, 143, 131 -
#' 143, 2013.
#'
#' T. W. Yee and A. G. Stephenson, Vector generalized linear and additive
#' extreme value models, Extremes, 10, 1 -- 19, 2007.
#' @keywords models
#' @examples
#'
#'   \donttest{
#'   #mod <- evm(rain, th=30)
#'   #mod
#'   #par(mfrow=c(2, 2))
#'   #plot(mod)
#'   }
#'
#'   \donttest{
#'   mod <- evm(rain, th=30, method="sim")
#'   par(mfrow=c(3, 2))
#'   plot(mod)
#'   }
#'
#'   \donttest{
#'   mod <- evm(SeaLevel, data=portpirie, family=gev)
#'   mod
#'   plot(mod)
#'   }
#'
#'   \donttest{
#'   mod <- evm(SeaLevel, data=portpirie, family=gev, method="sim")
#'   par(mfrow=c(3, 3))
#'   plot(mod)
#'   }
#'
#' @export evm
evm <- function(y, data, family=gpd, ...){
  theCall <- match.call()
  if (!missing(data)) {
      y <- ifelse(deparse(substitute(y))== "substitute(y)", deparse(y),deparse(substitute(y)))
      y <- formula(paste(y, "~ 1"))
      y <- model.response(model.frame(y, data=data))
  }
  UseMethod("evm", y)
}

#' @rdname evm
#' @export
evm.default <-
function (y, data, family=gpd, th= -Inf, qu,
          ..., # arguments specific to family such as phi = ~ 1
          penalty = NULL, prior = "gaussian",
          method = "optimize", cov="observed",
          start = NULL, priorParameters = NULL,
          maxit = 10000, trace=NULL,
          iter = 40500, burn=500, thin = 4, chains = 1,
          proposal.dist = c("gaussian", "cauchy"),
          jump.cov, jump.const=NULL,
          R=1000, cores=NULL, export=NULL, verbose=TRUE) {

    modelParameters <- texmexParameters(theCall, family,...)

    ##################### Sort out method, penalty/prior, trace...

    method <- texmexMethod(method)
    prior <- texmexPrior(prior, penalty, method, priorParameters)

    trace <- texmexTrace(trace, method)
    otrace <- trace[1]; trace <- trace[2]

    ############################## Construct data to use...

    if (missing(data)){ data <- NULL }
    else { y <- deparse(substitute(y)) }

    # Get list containing response (y) and design matrix for each parameter
    modelData <- texmexPrepareData(y, data, modelParameters)

    if (missing(th) & !missing(qu)) {
        th <- quantile(modelData$y, qu)
    }

    if (!is.finite(th)){ rate <- 1 }
    else { rate <- mean(modelData$y > th) }

    modelData <- texmexThresholdData(th, modelData)

    ###################### If family does not give info matrix...
    if (is.null(family$info)){ cov <- "numeric" }

    ###################### Check and sort out prior parameters...
    priorParameters <- texmexPriorParameters(prior, priorParameters, modelData)

    ################################## Do the optimization....
    o <- evmFit(data = modelData, family=family, ..., th=th,
                 prior=prior,
                 start=start, hessian = cov == "numeric",
                 priorParameters = priorParameters,
                 maxit = maxit, trace = otrace)

    if (o$convergence != 0 | o$value == 10^6) {
        warning("Non-convergence in evm.default")
    }

    ##### Construct object containing the penalized likelihood estimates
    o <- constructEVM(o, family, ..., th=th, rate=rate, prior=prior,
                      modelParameters=modelParameters, call=theCall,
                      modelData=modelData, data=data,
                      priorParameters=priorParameters, cov=cov)

    #### Simulate from posteriors....
    if (method == "s"){
        proposal.dist <- match.arg(proposal.dist)
        o <- evmSim(o, priorParameters=priorParameters,
                    prop.dist=proposal.dist,
                    jump.const=jump.const, jump.cov=jump.cov,
                    iter=iter, start=start, verbose=verbose,
                    thin=thin, burn=burn, chains=chains,
                    export=export,
                    trace=trace, theCall)
    } # Close else
    else if (method == "b"){
        o <- evmBoot(o, R=R, cores=cores, export=export)
    }

    o
}

