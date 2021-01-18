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

#' @title Method parameters for generating imptree objects
#' 
#' @description Initializing and validating the essential probability method specific parameters
#' 
#' @param args Named list containing the arguments to be processed.
#' May be \code{NULL} for default values. See details.
#' @param method Probability method as character, as supplied to \code{\link{imptree}}.
#' 
#' @return A list containing the sanitized and validated parameters.
#'
#' @details
#' \code{imptree_params()} is not exported into the user's namespace.
#' 
#' For all methods \code{args} takes the following inputs:
#' \itemize{
#' \item{s}{: Hyperparamter of the imprecise Dirichlet model
#'   (\code{s >= 0}), see below.}
#' \item{correction}{: Entropy correction to be carried out
#'   (Default \code{"no"}), see below.}
#' \item{splitmetric}{: Split criterion to use
#'   (Default \code{"globalmax"}), see below.}
#' }
#' The hyperparamter \code{s} of the imprecise Dirchlet model (IDM) may
#' be given as any non-negative value.  It defines the impression the locally
#' applied IDMs introduce. With increasing values of \code{s} more impression is
#' added. For \code{s=0} the IDM collapses to a precise Dirichlet model.
#' This value is ignored for \code{method = "NPI"}.
#' 
#' To account for a varying number of categories of the splitting candidates
#' Strobl proposed the use of a correction based on the Miller-entropy
#' correction: \code{correction = "strobl"}.
#' In their work Abellan and Moral favoured for the IDM the use of a
#' generalized Hartley measure such that the final measure may be viewed as
#' measure of total uncertainty: \code{correction = "abellan"}.
#' This correction method is not available for \code{method = "NPI"}.
#' 
#' When deciding for split canditates a split criterion is applied.
#' \code{"globalmax"} splits on maximal entropy of local models (with a
#' global IDM parameter \code{s}).
#' For \code{"range"} the splitting variable is found by taking the whole
#' entropy interval into account.
#' \code{localmax} is only available for IDM and split on maximal entropy,
#' however with \code{s} dependent on the number of missing values in the class
#' variable in the node 
#'  
#' @author Paul Fink \email{Paul.Fink@@stat.uni-muenchen.de}
#' 
#' @seealso \code{\link{imptree}}, \code{\link{imptree_control}}
#' 
#' @examples
#' ## Note: 
#' ## The function is used internally by imptree (not exported).
#' 
#' ## default constructed for method IDM
#' imptree:::imptree_params(NULL, method = "IDM")
#' 
#' ## passing arguments as list ('s' is not required for 'NPI')
#' imptree:::imptree_params(args = list(correction = "strobl", 
#'                                      splitmetric = "globalmax"),
#'                          method = "NPI")
#'
#' @keywords internal
imptree_params <- function(args, method) {

  # we expect a list or missing or null for 'args'
  if(missing(args) || is.null(args)) {
    args <- list(s = 1, correction = 0L, splitmetric = 0L)
  } else if (is.list(args)) {
    if (is.null(names(args))) {
      stop("the list supplied to 'method.param' must have names",
           domain ="R-imptree")
    }
    temp <- pmatch(names(args), c("s", "correction", "splitmetric"),
                   nomatch = 0L)
    if (any(temp == 0L)) {
      stop(gettextf("'method.param' component not matched: %s",
        names(args)[temp == 0L], domain ="R-imptree"))
    }
    names(args) <- c("s", "correction", "splitmetric")[temp]
    
    # check for IDM s
    if (is.null(args[["s"]])) {
      args$s <- 1
    } else if (method == "IDM" && args[["s"]] <= 0) {
      stop(gettextf("value of 's' (%.3f) must be strictly positive",
                    args$s, domain ="R-imptree"))
    }

    # check for 'correction'
    if (is.null(args[["correction"]])) {
      args$correction <- 0L
    } else {
      choices <- switch(substr(method, 1, 3), # remove the 'approx'
                        IDM = c("no", "strobl", "abellan"),
                        NPI = c("no", "strobl"))
      temp <- args[["correction"]][1]
      i <- pmatch(temp, choices, nomatch = 0L)
      if (i == 0L) {
        stop(gettextf("'correction' should be one of %s",
                      paste(dQuote(choices), collapse = ", "),
                      domain ="R-imptree"))
      }
      args$correction <- as.integer(i - 1)
    }
    # check for 'splitmetric'
    if (is.null(args[["splitmetric"]])) {
      args$splitmetric <- 0L
    } else {
      choices <- switch(substr(method, 1, 3), # remove the 'approx'
                        IDM = c("globalmax", "range"),# "localmax"),
                        NPI = c("globalmax", "range"))
      temp <- args[["splitmetric"]][1]
      i <- pmatch(temp, choices, nomatch = 0L)
      if (i == 0L) {
        stop(gettextf("'splitmetric' should be one of %s",
                      paste(dQuote(choices), collapse = ", "),
                      domain ="R-imptree"))
      }
      args$splitmetric <- as.integer(i - 1)
    }
  } else {
    stop("the argument supplied to 'method.param' must be a list or NULL for default values",
         domain ="R-imptree")
  }
  args
}
