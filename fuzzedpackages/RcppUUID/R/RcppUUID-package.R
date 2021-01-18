#' @title Generating Universally Unique Identificators
#'
#' @description
#' Provides functions to generating a vector of Universally Unique Identifiers (UUID).
#' Used implementation from the Boost C++ library. Supported random (version 4) and name (version 5)
#' UUIDs. UUIDs generation are parallelized by OpenMP.
#'
#' @name RcppUUID
#' @docType package
#'
#' @useDynLib RcppUUID, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
"_PACKAGE"
