#' Summary method for \code{pandemicEstimated} objects
#'
#' The summary method for \code{pandemicEstimated} object of class S3 displays a compact summary of the
#' fitted model. See the \strong{Details} section below for descriptions of the different components of the printed
#' output.
#'
#'
#' @method summary pandemicEstimated
#' @templateVar pandemicEstimatedArg object
#' @param object An object of S3 class \code{\link{pandemicEstimated-objects}}.
#' @param probs A numeric vector of quantiles of interest. The default is
#' \code{c(0.025,0.5,0.975)}.
#' @param digits Number of digits to use for formatting numbers.
#' @param info TRUE or FALSE: more details for output interpretation. The Default is TRUE.
#' @param ... Currently unused.
#' @details
#'
#' \subsection{Point estimates}{
#' Mean and Quantiles computed from simulations.
#' For models fit using MCMC ("sampling", this is default algorithim of \code{pandemic_model} function),  the posterior sample
#'  is used. For others estimation algorithm see \code{\link[rstan]{sampling}}  (\pkg{rstan} package).
#' }
#'
#' \subsection{Convergence and efficiency diagnostics for Markov Chains}{
#'
#' Included in the summary are: split effective sample sizes (n_eff), Monte Carlo standard errors
#' (se_mean) and split Rhats.
#'
#'
#' The Monte Carlo standard error provides relevant information for a posterior sample with more than one chain.
#'
#'
#' The R-hat convergence diagnostic compares the
#' between- and within-chain estimates for model parameters and other univariate
#' quantities of interest. If chains have not mixed well (ie, the between- and
#' within-chain estimates don't agree), R-hat is larger than 1.
#' We recommend running at least four chains by default and only using the
#' sample if R-hat is less than 1.05.
#'
#' }
#' \subsection{covidLPconfig}{
#' This subsection shows the main input settings used by the fitted model, and indicates whether default settings
#' of the CovidLP app (\url{http://est.ufmg.br/covidlp/home/en/})
#' were used (\code{covidLPconfig = TRUE} or \code{FALSE}).
#' Check the default settings of the CovidLP app in \code{\link{pandemic_model}}.
#' }
#' \subsection{Priors}{
#'
#' A list with information about the prior distributions used and model restrictions (if there are any).
#' For more information go to \code{\link{models}}.
#'
#'
#' }
#'
#' @return The summary method returns an object of class "summary.pandemicEstimated",
#'  which is a list with two arrays of summary statistics and diagnostics and
#'  others informations for use by the print method. The print method for summary.pandemicEstimated objects is called for its side effect
#'  and just returns its input.
#'
#'
#'@examples
#'\dontrun{
#' Y0=load_covid(country_name="Brazil",state_name="SP",last_date='2020-04-25')
#' output0=pandemic_model(Y0)
#' s=summary(output0)
#' s          #print method for summary.pandemicEstimated
#' names(s)  #see output list elements of summary
#' }
#'
#' @importMethodsFrom rstan summary
#'
#' @export
summary.pandemicEstimated=function(object,probs=c(0.025,0.5,0.975),digits=3,info=TRUE,...){

  if(class(digits)!= "numeric" | class(info)!= "logical"){
    stop("error in 'digits' or 'info'. View ?summary.pandemicEstimated")
  }

  # display parameters selection:
  pars=names(object$fit)[which( !grepl("mu",names(object$fit)) & !grepl("lp__",names(object$fit)) )]

  tab=rstan::summary(object$fit,pars=pars,probs=probs)$summary

  tab1=rstan::summary(object$fit,pars="mu",probs=probs)$summary
  takeout=which(colnames(tab1)=="n_eff" | colnames(tab1)=="Rhat" | colnames(tab1)=="se_mean") #convergence diagnosis does not make sense for mu
  tab1=tab1[,-takeout] #takeout convergence statistics for mu

  out=list(
    tab.parameters=tab,
  tab.fitted_values=tab1,
  cases.type=object$cases.type,
  seasonal_effect=object$seasonal_effect,
  n_waves=object$n_waves,
  location=object$Y$name,
  print.digits=digits,
  info=info,
  n=nrow(object$Y$data),
  covidLPconfig=object$config.inputs$covidLPconfig,
  init1=object$config.inputs$use_inputs$init[[1]],
  warmup=object$config.inputs$use_inputs$warmup,
  thin=object$config.inputs$use_inputs$thin,
  sample_size=object$config.inputs$use_inputs$sample_size,
  number_chains=object$config.inputs$use_inputs$number_chains,
  p=object$config.inputs$use_inputs$p
  )

  class(out)="summary.pandemicEstimated"

  return(out)

}






