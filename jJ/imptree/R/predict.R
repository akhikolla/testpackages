# Copyright (C) 2012 - 2018  Paul Fink
#
# This file is part of imptree.
#
# imptree is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# imptree is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with imptree.  If not, see <https://www.gnu.org/licenses/>.

#' @title Classification with Imprecise Probabilities
#' 
#' @description Prediction of \code{imptree} objects
#' 
#' @param object An object of class \code{imptree}. See details.
#' @param data Data.frame containing observations to be predicted.
#' If \code{NULL} the observations in the training set of \code{"object"}
#' are employed. 
#' @param dominance Dominace criterion to be applied when predicting
#' classes. This may either be \code{"strong"} (default) or \code{"max"}.
#' See details.
#' @param utility Utility for the utility based accuracy measure for a 
#' vacuous prediction result (default: 0.65).
#' @param \dots Additional arguments for data. May be \code{"weights"},
#' \code{"subset"}, \code{"na.action"}, any further are discarded.
#' 
#' @return \code{predict.imptree()} return an object of class 
#' \code{evaluation_imptree}, which is a named list containing 
#' predicted classes, predicted probability distribution and  accuracy 
#' evaluation
#' \item{probintlist}{List of the imprecise probability distributions of the
#' class variable. One matrix per observation in the test data.}
#' \item{classes}{Predicted class(es) of the observations as boolean matrix}
#' \item{evaluation}{Result of accuracy evaluation
#' \itemize{
#' \item{nObs}{: Number of observations}
#' \item{deter}{: Determinacy}
#' \item{nObsIndet}{: Number of observations with indeterminate prediction}
#' \item{indetSize}{: Average number of classes when predicting 
#' indeterminate (\code{NA} when no indeterminate observation)}
#' \item{acc_single}{: Single-set accuracy (\code{NA} when no determinate 
#' observation)}
#' \item{acc_set}{: Set-accuracy (\code{NA} when no indeterminate observation)}
#' \item{acc_disc}{: Discounted-accuracy}
#' \item{acc_util}{: Utility based (discounted) accuracy}
#' }}
#' 
#' @details
#' This function carries out the prediction of an imprecise tree. 
#' An existence check on the stored C++ object reference is carried out 
#' at first. If the reference is not valid the original call
#' for \code{"object"} is printed as error.
#' 
#' There are currently 2 different dominance criteria available:
#' \describe{
#' \item{max}{Maximum frequency criterion. Dominance is decided only
#' by the upper bound of the probability interval, ie. a state \eqn{C_i} is 
#' dominated if there exists any \eqn{j \neq i}{j != i} with
#' \eqn{u(C_i) < u(C_j)}}
#' \item{strong}{Interval dominance criterion. For the IDM it
#' coincides with the strong dominance criterion. Here a state
#' \eqn{C_i} is dominated if there exists any \eqn{j \neq i}{j != i} 
#' with \eqn{u(C_i) < l(C_j)}}
#' }
#' @author Paul Fink \email{Paul.Fink@@stat.uni-muenchen.de}
#' 
#' @seealso \code{\link{imptree}}, \code{\link{node_imptree}}
#' 
#' @keywords tree
#'
#' @examples
#' data("carEvaluation")
#' 
#' ## create a tree with IDM (s=1) to full size
#' ## carEvaluation, leaving the first 10 observations out
#' ip <- imptree(acceptance~., data = carEvaluation[-(1:10),], 
#'   method="IDM", method.param = list(splitmetric = "globalmax", s = 1), 
#'   control = list(depth = NULL, minbucket = 1))
#' 
#' ## predict the first 10 observations with 'max' dominance
#' pp <- predict(ip, dominance = "max", data = carEvaluation[(1:10),])
#' print(pp)
#' pp$classes                ## predicted classes as logical matrix
#' 
#' ## predict the first 10 observations with 'strong' dominance and
#' ## use a different level of utility
#' predict(ip, dominance = "strong", data = carEvaluation[(1:10),],
#'         utility = 0.5)
#'
#' @export
predict.imptree <- function(object, data, dominance = c("strong", "max"),
                            utility = 0.65, ...) {

  # Are the C++ object references still stored in the R object?
  tryCatch({hasRoot_cpp(object$tree)},  error = function(e) {
    stop(gettextf("reference to tree is not valid; see element \"call\" of '%s' for recreation", 
                 deparse(substitute(object)), domain = "R-imptree"))
  })
  
  # checking for valid utility input
  if(0 >= utility || 1 <= utility) {
    stop(gettextf("value of 'utility' (%f) needs to be in (0,1)", 
                  utility, domain = "R-imptree"))
  }
  
  # matching dominance type to C++ as int 
  dominance <- match.arg(dominance)
  domint <- pmatch(dominance, c("strong", "max")) - 1
  
  predcontrol <- list(utility = as.double(utility),
                      dominance = as.integer(domint))
  
  # In case of missing data, evaluation and predicting on training data
  if(missing(data)) {
    testdata <- object$traindata
  } else {
    testdata <- prepare_data(object, data = data, ...)
    
  }
  evaluation <- predict_cpp(object$tree, testdata, predcontrol)

  # Transposing the class matrix so observations are in rows
  evaluation$classes <- t(evaluation$classes)  

  # Setting the classlabels on the prediction matrix and the Probinterval list
  classlabels  <- attr(object$traindata, "labels")[[1]]

  colnames(evaluation$classes) <- classlabels
	
  evaluation$probintlist <- lapply(evaluation$probintlist, 
                                   function(o) {
    colnames(o) <- classlabels
    o
  })

  # adding the call parameters
  attr(evaluation, "dominance") <- dominance
  attr(evaluation, "utility") <- utility
  class(evaluation) <- c("evaluation_imptree", class(evaluation))
  evaluation
}

#' @rdname predict.imptree
#' @param x an object of class \code{evaluation_imptree}
#' @return The printing function returns the
#' \code{evaluation_imptree} object invisibly.
#' @export
print.evaluation_imptree <- function(x, ...) {
  
  meval <- matrix(unlist(x$evaluation)[-c(1,3)], ncol = 1)
  dimnames(meval) <- list(gettextf(c("Determinacy",
                                     "Average indeterminate size",
                                     "Single-Set Accuracy",
                                     "Set-Accuracy",
                                     "Discounted Accuracy",
                                     "%.2f utility based Accuracy"),
                                   attr(x, "utility"), 
                                   domain ="R-imptree"),
  "")
  cat(gettextf(
    "Accuracy achieved on %d observations in testing data, based on '%s' dominance:\n",
              x$evaluation$nObs,
              attr(x, "dominance"), domain = "R-imptree"))
  if((detObs <- x$evaluation$nObs - x$evaluation$nObsIndet) > 0) {
    cat(gettextf("\t%d determinate predictions\n", 
                 as.integer(detObs), 
                 domain ="R-imptree"))
  }
  if(x$evaluation$nObsIndet > 0) {
    cat(gettextf("\t%d indeterminate predictions\n",
                 as.integer(x$evaluation$nObsIndet), 
                 domain ="R-imptree"))
  }
  print(meval)
  invisible(x)
}
