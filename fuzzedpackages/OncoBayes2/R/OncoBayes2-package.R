#' OncoBayes2
#'
#' Bayesian logistic regression model with optional
#' EXchangeability-NonEXchangeability parameter modelling for flexible
#' borrowing from historical or concurrent data-sources. The safety
#' model can guide dose-escalation decisions for adaptive Oncology
#' phase I dose-escalation trials which involve an arbitrary number of
#' drugs.
#'
#' @section Global Options:
#'
#' \tabular{lcl}{
#' Option \tab Default \tab Description \cr
#' \code{OncoBayes2.MC.warmup} \tab 1000 \tab MCMC warmup iterations \cr
#' \code{OncoBayes2.MC.iter} \tab 2000 \tab total MCMC iterations \cr
#' \code{OncoBayes2.MC.chains} \tab 4 \tab MCMC chains\cr
#' \code{OncoBayes2.MC.thin} \tab 1 \tab MCMC thinning \cr
#' \code{OncoBayes2.MC.control} \tab \code{list(adapt_delta=0.99,} \tab sets \code{control} argument for Stan call\cr
#'  \tab \code{stepsize=0.1}) \tab \cr
#' \code{OncoBayes2.abbreviate.min} \tab 0 \tab Minimal length of variable names \cr
#'    \tab \tab when abbreviating variable names. \cr
#'    \tab \tab The default 0 disables abbreviation.\cr
#' }
#'
#' @template ref-mac
#' @template ref-exnex
#' @template ref-critical_aspects
#' @template ref-bayesindustry
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
#' @name OncoBayes2
#' @alias OncoBayes2
#' @docType package
#' @useDynLib OncoBayes2, .registration = TRUE
#' @import methods
#' @importFrom rstan sampling get_sampler_params summary stanc_builder
#' @importFrom utils capture.output modifyList combn
#' @importFrom stats delete.response ftable median model.frame model.matrix model.response quantile rbinom sd terms model.matrix.default setNames update update.default .getXlevels as.formula na.fail qlogis
#' @import methods
#' @import assertthat
#' @import checkmate
#' @import Formula
#' @import Rcpp
#' @import rstantools
#' @import dplyr
#' @import tidyr
#' @importFrom tibble as_tibble
#' @import abind
#' @export posterior_linpred posterior_predict posterior_interval
#' @export predictive_interval prior_summary nsamples
#'
#'
NULL
