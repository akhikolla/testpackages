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

#' @title Classification Trees with Imprecise Probabilities
#' 
#' @name imptree
#' 
#' @description \code{imptree} implements Abellan and Moral's tree 
#' algorithm (based on Quinlans ID3) for classification. It
#' employes either the imprecise Dirichlet model (IDM) or
#' nonparametric predictive inference (NPI) to generate the
#' imprecise probability distribution of the classification variable
#' within a node.
#' 
#' @param formula Formula describing the strucutre
#' (class variable ~ featutre variables).
#' Any interaction terms trigger an error.
#' @param data Data.frame to evaluate supplied formula on.
#' If not provided the the formula is evaluated 
#' on the calling environment
#' @param weights Individual weight of the observations
#' (default: 1 to each).
#' \emph{This argument is ignored at the moment.}
#' @param control A named (partial) list according to the result of
#' \code{\link{imptree_control}}.
#' @param method Method applied for calculating the probability
#' intervals of the class probability. \code{"IDM"} for the imprecise
#' Dirichlet model (default), \code{"NPI"} for use of the 
#' nonparametric predictive inference approach and \code{"NPIapprox"}
#' for use of the approximate algorithm obtaining maximal entropy of
#' NPI generated probability intervals.
#' @param method.param Named list providing the method specific 
#' parameters. See \code{\link{imptree_params}}.
#' @param \dots optional parameters to be passed to the main function
#' \code{imptree.formula} or to the call of
#' \code{\link{imptree_control}}.
#'
# This part is not yet implemented. Will be done hopefully in a later version
# @details
# A multi-state observation of a variable when treated epistemically, all 
# the states are possible values of the observations, however in case of
# ontological interpretation each unique multi-state is assumed as a further
# single state. The a multi-state observation is characterized by a factor
# with the levels for a multi-state consisting of the unique single stated
# separated by a colon \code{":"}.
#  
#' @return An object of class \code{imptree}, which is a list
#' with the following components:
#' \item{call}{Original call to \code{imptree}}
#' \item{tree}{Object reference to the underlying C++ tree object.}
#' \item{train}{Training data in the form required by the 
#' workhorse C++ function.\cr
#' It is an integer matrix containing the internal factor
#' representations, adjusted for the C++ specific indexing
#' starting at 0 and not at 1 as in R.
#' Further attributes of the matrix, hold the names of the variables,
#' the C++ adjusted index of the classification variabe, as well as
#' the levels and number of levels for each variable.}
#' \item{formula}{The formula describing the data structure}
#' 
#' @references Abell\ifelse{latex}{\out{\'{a}}}{\ifelse{html}{\out{&aacute;}}{a}}n,
#' J. and Moral, S. (2005), Upper entropy of credal sets. Applications to 
#' credal classification, \emph{International Journal of Approximate Reasoning}
#' \bold{39}, 235--255.
#' @references Strobl, C. (2005), Variable Selection in Classification Trees Based on
#' Imprecise Probabilities, \emph{ISIPTA'05: Proceedings of the Fourth
#' International Symposium on Imprecise Probabilities and Their Applications},
#' 339--348.
#' @references Baker, R. M. (2010), \emph{Multinomial Nonparametric Predictive Inference:
#' Selection, Classification and Subcategory Data}.
#'  
#' @author Paul Fink \email{Paul.Fink@@stat.uni-muenchen.de},
#' based on algorithms by 
#' J. Abell\ifelse{latex}{\out{\'{a}}}{\ifelse{html}{\out{&aacute;}}{a}}n
#' and S. Moral for the IDM and R. M. Baker for the NPI approach.
#' 
#' @seealso \code{\link{predict.imptree}} for prediction,
#' \code{\link{summary.imptree}} for summary information, 
#' \code{\link{imptree_params}} and \code{\link{imptree_control}} for
#' arguments controlling the creation, \code{\link{node_imptree}} for
#' accessing a specific node in the tree
#' 
#' @keywords tree
#'
#' @examples
#' data("carEvaluation")
#' 
#' ## create a tree with IDM (s=1) to full size on
#' ## carEvaluation, leaving the first 10 observations out
#' imptree(acceptance~., data = carEvaluation[-(1:10),], 
#'   method="IDM", method.param = list(splitmetric = "globalmax", s = 1), 
#'   control = list(depth = NULL, minbucket = 1)) # control args as list
#' 
#' ## same setting as above, now passing control args in '...'
#' imptree(acceptance~., data = carEvaluation[-(1:10),], 
#'   method="IDM", method.param = list(splitmetric = "globalmax", s = 1), 
#'   depth = NULL, minbucket = 1)
#' 
#' @importFrom stats terms
#' @export
imptree.formula <- function(formula, data = NULL,
                            weights, control, method = c("IDM", "NPI", "NPIapprox"),
                            method.param, ...) {

  # making sure a valid formula is supplied  
  if(missing(formula) || !inherits(formula, "formula")) {
    stop("'formula' is missing but mandatory", domain ="R-imptree")
  }
  
  # function call
  Call <- match.call()
  
  # generating a suitable data set representation
  dataset <- prepare_data(object = formula, data = data,
                            weights = weights, ...)

  # checking for the calculation method
  method <- match.arg(method)
  methodint <- as.integer(pmatch(method, c("IDM", "NPI", "NPIapprox")) - 1)
  
  # checking for the method specific parameters
  method.param <- imptree_params(method.param, method)
  
  # initializing the control arguements such as minbucket, depth,
  # gamma and tbase, they may also be supplied via '...'
  if (missing(control)) {
    control <- NULL
  }
  controls <- imptree_control(splitmetric = method.param$splitmetric,
                              controlList = control, ...)

  # combine all controls and method parameters into a single list
  controlslist <- c(controls, method.param)
  controlslist$iptype <- methodint
  
  # init java tree 
  tree_cpp_object <- treebuilder_cpp(dataset, controlslist)

  # changing back the call to 'imptree'
  Call[[1]] <- as.name("imptree")
  # construction of the structure of object 'imptree'
  res <- list(call = Call,
              tree = tree_cpp_object,
              traindata = dataset,
              formula = formula(terms(formula, data = dataset)))
  # assigning class
  class(res) <- c("imptree")
  # return value
  res
}


#' @rdname imptree
#' @param x A data.frame or a matrix of feature variables.
#' The columns are required to be named.
#' @param y The classification variable as a factor.
#' @importFrom stats as.formula
#' @export
imptree.default <- function(x, y, ...) {
  
  if(missing(y) || missing(x)) {
    stop("arguments 'y' and 'x' are mandatory", 
         domain ="R-imptree")
  }
  
  if(is.null(nam <- colnames(x))) {
    stop(gettextf("'%s' must contain column names", 
                 deparse(substitute(x)), domain ="R-imptree"))
  }
  formula <- as.formula(paste0(deparse(substitute(y)),"~",
                               paste(nam, collapse = "+")))
  res <- imptree.formula(formula,
                         data = data.frame(deparse(substitute(y)) , x),
                         ...)
  call <- match.call()
  call[[1]] <- as.name("imptree")
  res$Call <- call
  res  
}


#' @rdname imptree
#' @export
imptree <- function(x, ...) {
  UseMethod("imptree")
}
