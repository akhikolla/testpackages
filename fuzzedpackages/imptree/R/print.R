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
#' @description Printing the \code{imptree} object to console
#' 
#' @param x Object of class \code{imptree}. See details.
#' @param digits a non-null value for digits specifies the minimum number
#' of significant digits to be printed in values. The default uses 
#' \code{\link[base]{getOption}("digits")}. Non-integer values will be rounded down,
#' and only values greater than or equal to 1 and 
#' no greater than 17 are accepted.
#' @param sep Separator between the displayed IPDistribution objects.
#' (Default: \code{'\t'})
#' @param \dots Additional arguments; ignored at the moment
#' 
#' @return Returns the calling object invisible.
#' 
#' @details
#' An existence check on the stored C++ object reference is carried out 
#' at first. If the reference is not valid the original call
#' for \code{"object"} is printed as error.
#' 
#' For a more detailed summary of the tree \code{\link{summary.imptree}}.
#'
#' @author Paul Fink \email{Paul.Fink@@stat.uni-muenchen.de}
#' 
#' @seealso \code{\link{imptree}}, \code{\link{summary.imptree}}
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
#' ip                        ## standard printing; same as 'print(ip)'
#' print(ip, sep = ";")      ## probability intervals are separated by ';'
#' 
#' @export
print.imptree <- function(x, digits = getOption("digits"), sep = "\t", ...){
  
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  # Are the C++ object references still stored in the R object?
  tryCatch({hasRoot_cpp(x$tree)},  error = function(e) {
    stop(gettextf("reference to tree is not valid; see element \"call\" of '%s' for recreation", 
                 deparse(substitute(x)), domain = "R-imptree"))
  })
  digits <- floor(digits);
  if(digits < 1 || digits > 17) {
    stop("invalid 'digits' argument", domain = "R-imptree")
  }
  
  treePrint_cpp(x$tree, nsmall = digits, sep = sep)
  invisible(x)
}
