# Copyright (C) 2018  Paul Fink
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

#' @title Various method around IPIntervals
#' 
#' @description Calculation of probability intervals, 
#' and their maximal and minimal entropy
#' 
#' @param table integer vector of absolute frequencies
#' @param iptype method for calculating the probability
#' intervals of \code{table}. \code{"IDM"} for the imprecise
#' Dirichlet model (default), \code{"NPI"} for use of the 
#' nonparametric predictive inference approach and \code{"NPIapprox"}
#' for use of the approximate algorithm obtaining maximal entropy of
#' NPI generated probability intervals.
#' @param entropymin Calculation of one distribution with minimal 
#' entropy, including the actual value of the minimal entropy
#' (default: \code{TRUE})
#' @param entropymax Calculation of the distribution with maximal 
#' entropy, including the actual value of the maximal entropy
#' (default: \code{TRUE})
#' @param correction Entropy correction to be carried out,
#' ignorned if \code{(entropymin || entropymax) == FALSE}
#' (default \code{"no"}), see \code{\link{imptree_params}}
#' @param s Hyperparamter of the IDM (\code{s >= 0}),
#' see \code{\link{imptree_params}} 
#' (ignored for \code{iptype == "NPI"})
#' 
#' @return A list with 5 named entries:
#' \item{probint}{matrix with 3 rows and \code{length(table)}
#' columns: in the rows are the abosulte frequencies, the lower 
#' bound (\code{"lower"}) and the upper bound (\code{"upper"}) 
#' of the event-wise probabilities.}
#' \item{maxEntDist}{The (unique) probability distribution with
#' maximal entropy}
#' \item{maxEntCorr}{The value of the (corrected) maximal entropy}
#' \item{minEntDist}{A probability distribution with minimal
#' entropy, as it is not necessarily unqiue there may be others}
#' \item{minEntCorr}{The value of the (corrected) minimal entropy}
#' 
#' @author Paul Fink \email{Paul.Fink@@stat.uni-muenchen.de}
#' 
#' @seealso \code{\link{imptree_params}}
#' 
#' @examples
#' ## Artificial vector of absolute frequencies
#' obs <- c(a = 1,b = 2, c = 10, d = 30, e = 5)
#'
#' ## probability interval by NPI, including only information on the
#' ## mininum entropy distribution, using no entropy correction
#' probInterval(obs, iptype = "NPI", entropymax = FALSE)
#'
#' ## probability interval by IDM, including information on the
#' ## minimum and maximum entropy distribution with s = 2 and correction
#' ## according to 'strobl'
#' probInterval(obs, iptype = "IDM", correction = "strobl", s = 2)
#' 
#' @export
probInterval <- function(table, 
                         iptype = c("IDM", "NPI", "NPIapprox"), 
                         entropymin = TRUE, entropymax = TRUE,
                         correction = c("no", "strobl", "abellan"),
                         s = 1) {
  choices <- c("IDM", "NPI", "NPIapprox")
  ipt <- pmatch(iptype[1], choices, nomatch = 0L)
  if(ipt == 0L) {
    stop(sprintf(gettext("'iptype' should be one of %s"), 
                 paste(dQuote(choices), collapse = ", ")))
  }
  ipt <- as.integer(ipt - 1)
  
  if (ipt == 0L && s <= 0) {
    stop(sprintf(
      gettext("value of 's' (%f) must be strictly positive"), s))
  }
  if(entropymin || entropymax) {
    choices <- if(ipt == 0L) {
      c("no", "strobl", "abellan")
    } else {
      c("no", "strobl")
    }
    entcorr <- pmatch(correction[1], choices, nomatch = 0L)
    if (entcorr == 0L) {
      stop(sprintf(gettext("'correction' should be one of %s"),
                   paste(dQuote(choices), collapse = ", ")),
           domain = "imptree")
    }
    entcorr <- as.integer(entcorr - 1)
  } else {
    entcorr <- 0L
  }
  
  config <- list(iptype = ipt, s = s, correction = entcorr)
  createProbIntInformation_cpp(as.integer(table), config, 
                               as.logical(entropymin), 
                               as.logical(entropymax))
}
