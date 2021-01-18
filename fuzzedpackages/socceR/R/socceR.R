#' Evaluating sport tournament predictions
#' 
#' Functions for evaluating sport tournament predictions, the tournament rank probability score, and working with models for prediction sport matches.
#' 
#' @docType package
#' @author Claus Ekstrom <ekstrom@@sund.ku.dk> 
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib socceR
#' @name socceR
NULL  

#' FIFA 2018 prediction matrices
#'
#' A list containing five predictions for the FIFA 2018 World Cup.
#'
#' @format A list with 5 predictions (each a 7 by 32 matrix) containing the predictions probabilities of 1st, 2nd, 3rd, 4th, 5th-8th, 9th-12th, and 17th-32nd place.
#' \describe{
#'   \item{flat}{A prediction with equal probability of winning for all teams}
#'   \item{ekstrom1}{Ekstrom's prediction (based on the Skellam distribution)}
#'   \item{ekstrom2}{Ekstrom's prediction (based on the ELO rankings)}
#'   \item{GLSE1}{Prediction of Groll et all}
#'   \item{GLSE2}{Updated prediction of Groll et all}
#' }
#' @source \url{http://sandsynligvis.dk/2018/08/03/world-cup-prediction-winners/}
#'
"fifa2018"

#' FIFA 2018 end results
#'
#' A named vector sorted in the ranking of the teams in the FIFA 2018 World Cup. The value correspond to the corresponding columns in the prediction matrices of fifa2018 
#'
#' @format A vector of the final rankings
#' @source \url{http://sandsynligvis.dk/2018/08/03/world-cup-prediction-winners/}
#'
"fifa2018result"
