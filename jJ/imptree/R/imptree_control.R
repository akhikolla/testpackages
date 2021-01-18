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

#' @title Control parameters for generating
#' imptree objects
#' 
#' @description Initializing and validating 
#' the tree generation parameters
#' 
#' @param splitmetric Choosen split metric as integer:
#' \code{0} means \code{"globalmax"} and \code{1L} \code{"range"}, 
#' repectively. See \code{\link{imptree_params}}
#' @param controlList Named list containing the processed arguments.
#' See details.
#' @param tbase Value that needs to be at least attained to qualify
#' for splitting (default: 1)
#' @param gamma Weighting factor of the maximum entropy
#' (default: 1)
#' @param depth Integer limiting the tree to the given depth, with
#' \code{0} indicating to perform no splitting at all.
#' If not supplied, \code{NULL} (default) or negative the
#' tree is grown to maximal size, the latter triggering a warning.
#' @param minbucket Positive integer as minimal leaf size
#' (default: 1)
#' @param \dots Argument gobbling; is not processed
#' 
#' @return A list containing the options. Missing options are set
#' to their default value.
#' 
#' @details
#' The argument \code{controlList} may be a named list with names in
#' \code{c("tbase", "gamma", "depth", "minbucket")}
#' Any values in this list will overwrite those supplied in 
#' named arguments.
#' When \code{controlList = NULL} (default) only the supplied 
#' arguments are checked.
#' 
#' In case \code{controlList} contains an argument named
#' \code{splitmetric}, this will be ignored.
#' If \code{splitmetric} is \code{0L}, i.e. \code{"globalmax"}, 
#' the values for \code{gamma} and \code{tbase} are set to their 
#' default values, even if the user supplied different values.
#' 
#' @author Paul Fink \email{Paul.Fink@@stat.uni-muenchen.de}
#' 
#' @seealso \code{\link{imptree}}, \code{\link{imptree_params}}
#' 
#' @keywords tree
#' 
#' @examples
#' ## Check performed for splitmetric 'globalmax',
#' ## tbase' is default generated and 'gamma' is overwritten
#' ## (see Details), tree is grown to full depth and 
#' ## at least 5 observations are needed to be within each node
#' imptree_control(splitmetric = 0, gamma = 0.5,
#'                 depth = NULL, minbucket = 5)
#' 
#' ## Passing some control arguments in a list
#' ## As splitmetric is 'range', gamma is respected
#' imptree_control(splitmetric = 1, minbucket = 5,
#'                 controlList = list(gamma = 0.5, depth = NULL))
#'
#' @export
imptree_control <- function(splitmetric, controlList = NULL, tbase = 1, gamma = 1,
                            depth = NULL, minbucket = 1L, ...) {

  # save the method parameter list
  sm <- splitmetric
  
  # generation of the bare returned list
  clist <- list(depth = depth, minbucket = minbucket,
                tbase = tbase, gamma = gamma)
  
  # first we process add the agruments supplied in controlList
  if(!is.null(controlList)) {	
    clist[names(controlList)[names(controlList) %in% names(clist)]] <- 
      controlList[names(controlList) %in% names(clist)]
  }
  
  if(any(nal <- is.na(clist))) {
    stop(gettextf("no 'NA' permitted in %s", 
                 paste(sQuote(names(clist)[nal]), collapse = ", "),
                 domain ="R-imptree"))
  }

  # heading off now to test for inconsistencies
  if(sm) {

    # dealing inconsitencies in 'tbase'
    tbase <- as.double(clist[["tbase"]])
    if(tbase > 2 || tbase < -1) {
      stop(gettextf("value of 'tbase' (%.3f) must be between [-1,2]", 
                    tbase, domain ="R-imptree"))
    }
  
    # dealing inconsitencies in 'gamma'
    gamma <- as.double(clist[["gamma"]])
    if(gamma > 1 || gamma < 0) {
      stop(gettextf("value of 'gamma' (%.3f) must be in [0,1]",
                   gamma, domain ="R-imptree"))
    }
  } else {
  
    # if splitmetric is globalmax (i.e. 0L), then tbase and gamma are predefined
    clist$gamma <- 1
    clist$tbase <- 2
  }
  
  # dealing inconsitencies in 'depth'
  mydepth <- clist[["depth"]]
  if(!is.null(mydepth) && mydepth < 0L) {
    warning(gettextf("ignoring supplied 'depth'=%d and use default instead",
                    mydepth, domain ="R-imptree"))
    mydepth <- NULL
  }
  if(is.null(mydepth)) {
    mydepth <- as.integer(.Machine$integer.max)
  }
  clist$depth <- as.integer(mydepth)
  
  # dealing with inconsistency in 'minbucket'
  if((minbucket <- clist[["minbucket"]]) < 1L) {
    warning(gettextf("ignoring supplied 'minbucket'=%d and use default instead",
                    minbucket, domain ="R-imptree"))
    minbucket <- 1L
  }
  clist$minbucket <- as.integer(minbucket)
  clist
}
