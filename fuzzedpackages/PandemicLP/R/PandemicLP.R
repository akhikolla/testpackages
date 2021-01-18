#' PandemicLP: Modeling Pandemic Data
#'
#' \if{html}{\figure{logo.jpeg}{options: width="120px"}}
#'
#' The \pkg{PandemicLP} package provides five main functions that enable the user to make short and long-term
#' predictions of epidemic and pandemic count data.
#' \cr
#' This package is a result of the joint work by professors and graduate students from the Statistics
#' Department at Universidade Federal de Minas Gerais (UFMG). It originated as a challenge in a graduate course
#' after the suspension of classes due to Covid-19.
#' \cr
#' Theoretical foundation and more information about the project can be found in
#' \href{http://est.ufmg.br/covidlp/home/en/}{est.ufmg.br/covidlp/home/en}
#'
#' @author
#' Authors:
#' \itemize{
#' \item{Debora de Freitas Magalhaes}
#' \item{Marta Cristina Colozza Bianchi da Costa}
#' \item{Guido Alberti Moreira}
#' \item{Thais Pacheco Menezes}
#' \item{Marcos Oliveira Prates}
#' }
#'
#' @examples
#' \dontrun{
#' data = load_covid("Brazil")
#' plot(data)
#' estim = pandemic_model(data)
#' pred = posterior_predict(estim)
#' stats = pandemic_stats(pred)
#' plot(pred)}
#'
#' @docType package
#' @name PandemicLP
#' @useDynLib PandemicLP, .registration=TRUE
#' @import Rcpp
#' @import methods
#' @import dplyr
NULL
