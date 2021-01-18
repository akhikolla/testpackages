#' @title Simulation-based inference for infectious disease models
#'
#' @description Package implements various simulation-based inference routines for infectious
#' disease models.
#'
#' @details Provides some code to run simulations of state-space models, and then
#'  use these in the Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) 
#'  algorithm of Toni et al. (2009) <doi:10.1098/rsif.2008.0172> and a bootstrap particle
#'  filter based particle Markov chain Monte Carlo (PMCMC) algorithm 
#'  (Andrieu et al., 2010 <doi:10.1111/j.1467-9868.2009.00736.x>). 
#'  Also provides functions to plot and summarise the outputs.
#'
#' @docType package
#' @name SimBIID-package
#' @author Trevelyan J. McKinley <t.mckinley@@exeter.ac.uk>
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import ggplot2
#' @import tidyr
#' @import Rcpp 
#' @importFrom Rcpp evalCpp
#' @importFrom graphics plot
#' @importFrom stats window
#' @useDynLib SimBIID, .registration = TRUE
NULL

