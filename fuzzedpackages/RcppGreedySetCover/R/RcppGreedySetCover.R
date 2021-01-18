#' Fast Greedy Set Cover
#'
#' This package offers an implementation of the greedy algorithm for the set cover problem using Rcpp.
#' 
#' @references Vazirani, Vijay V. (2001), Approximation Algorithms, Springer
#' 
#' @section Set Cover Problem: 
#' We are given a universe of elements \emph{U} and \emph{F}, 
#' a collection of subsets from \emph{U}, covering \emph{U}; i.e. \emph{U} is 
#' in the union of the sets in \emph{F}.
#' The objective is to find \emph{A}, the smallest subcollection of \emph{F}, covering \emph{U}.
#' An important application is hospital placement, where the 
#' number of hospitals is minimized under the constraint that all residents are provided for.
#' 
#' The optimal solution to the problem is available via linear programming, 
#' however this is not a feasible solution for large 
#' problems due to the computational demands involved. 
#' A quick approximate solution is given by the greedy algorithm. 
#' The algorithm iterates the following steps until all elements are covered, 
#' starting from an empty \emph{A}:
#' \itemize{
#'   \item Add the largest set of uncovered elements to \emph{A}.
#'   \item Remove covered elements from \emph{F}.
#' }
#' This simple algorithm exhibits surprisingly good properties. 
#' For a nice introduction to the set cover problem and the greedy algorithm see Vazirani, 2001.
#' 
#' @docType package
#' @name RcppGreedySetCover-package
NULL