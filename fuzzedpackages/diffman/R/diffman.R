#' \code{diffman} package
#'
#' Package to manage differentiation
#' 
#' Differentiation is a technique to deduce information on small
#' agregates from a source of information that has been
#' disseminated according to two different nomenclatures. Knowing
#' the value of additive variables on agregates of observations
#' according to the nomenclatures, one can decude the value of
#' these variables on new agregates. If those agregates are 
#' smaller than a given threshold then we can considerer there
#' is a confidentiality break for the observations within those
#' agregates. It order to respect statistical secrecy for every
#' observations, one need to check for the possible problematic 
#' differentiations that a data user can compute. The objective
#' of this package is to provide a tool for detecting all
#' observations at risk of differentiation.
#'
#' @docType package
#' @name diffman
NULL

#' @useDynLib diffman
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#' @import data.table
NULL

## quiets concerns of R CMD check : the .'s that appear in data.table
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

