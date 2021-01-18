# Individual Based Models in R --------------------------------------------

#' @title Individual based models in R
#' @description Implementation of some (simple) Individual Based Models and methods 
#' to create new ones, particularly for population dynamics models (reproduction, 
#' mortality and movement). The basic operations for the simulations are 
#' implemented in Rcpp for speed.
#' @name ibm-package
#' @aliases ibm-package ibm
#' @docType package
#' @author Ricardo Oliveros-Ramos
#' <ricardo.oliveros@@gmail.com>
#' @keywords ibm, individual based models, agent based models
#' @examples
#' \dontrun{
#' set.seed(880820)
#' par = list(alpha=5e-4, beta=5e-4, r=0.1, m=0.05, D=list(N=8e-5, P=8e-5), 
#' L=list(N=0.2, P=0.2))
#' N0 = with(par, m/(2*beta*L$P))
#' P0 = with(par, r/(2*alpha*L$N))
#' par$initial = list(N=round(N0), P=round(P0))
#' sim = localLotkaVolterra(par, T=240, replicates=100, maxpop = 1e4)
#' plot(sim)
#' } 
NULL
