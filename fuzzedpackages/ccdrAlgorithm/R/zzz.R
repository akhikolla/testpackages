#
#  zzz.R
#  ccdrAlgorithm
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

### Imported from sparsebnUtils
# NOTE: These imports MUST go in this file -- they are ignored (why?) if they appear in other places
#' @importFrom sparsebnUtils reIndexC
#' @importFrom sparsebnUtils reIndexR
#' @importFrom sparsebnUtils get.adjacency.matrix
#' @importFrom sparsebnUtils num.nodes
#' @importFrom sparsebnUtils num.edges
#' @importFrom sparsebnUtils is.zero
#' @importFrom sparsebnUtils edgeList
#' @importFrom sparsebnUtils sparse
#' @importFrom sparsebnUtils is.edgeList
#' @importFrom sparsebnUtils to_igraph
#' @importFrom sparsebnUtils as.sparse
#' @importFrom sparsebnUtils is.sparse
#' @importFrom stats rnorm
#' @importFrom utils tail

.onAttach <- function(libname, pkgname){
    ### Only sparsebn needs a package startup message
    # packageStartupMessage("NOTE: This package is currently in a development state and may be unstable.\n Please report any bugs at https://github.com/itsrainingdata/ccdrAlgorithm/issues.")
}
