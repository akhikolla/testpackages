#
#  s3-SparseBlockMatrixR.R
#  ccdrAlgorithm
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2017 Bryon Aragam. All rights reserved.
#

#------------------------------------------------------------------------------#
# SparseBlockMatrixR S3 Class for R
#------------------------------------------------------------------------------#

#
# SparseBlockMatrixR S3 class skeleton
#
# Data
# * list rows
# * list vals
# * list blocks
# * numeric sigmas
# * integer start
#

#
# A convenience class to make sharing data between R and C++ easier. This class mimics the structure
#   of the C++ class 'SparseBlockMatrix' (note the name difference to differentiate the two) as a list in R,
#   which makes it easy to use Rcpp to pass a sparse structure between R and C++. This class is NOT intended
#   to be a general purpose sparse data structure for R; instead, it simply streamlines the connection between
#   R and C++.
#
# A SparseBlockMatrixR object consists of four main components:
#   1) rows - a list of integer vectors
#   2) vals - a list of numeric vectors
#   3) blocks - a list of integer vectors
#   4) sigmas - a numeric vector
#
# There is also a fifth component which identifies whether the indexing in the object begins at 0 or 1.
#   This is needed for bookkeeping and ensuring coherent translation between R and C++.
#   5) start - 0 or 1
#
# These components all exactly reflect their purpose in SparseBlockMatrix.h; see that file for more details.
#
#   NOTES:
#       1) Since C++ begins indexing at 0, we have included two functions for re-indexing between the two conventions,
#          defined in reIndexR and reIndexC. The 'start' flag ensures we keep track of where indexing begins.
#
#

### Need to import some generics from sparsebnUtils
#' @importFrom sparsebnUtils to_graphNEL sparse

#------------------------------------------------------------------------------#
# is.SparseBlockMatrixR
#
is.SparseBlockMatrixR <- function(x){
    inherits(x, "SparseBlockMatrixR")
} # END IS.SPARSEBLOCKMATRIXR

as.SparseBlockMatrixR <- function(x){
    SparseBlockMatrixR(x) # NOTE: S3 delegation is implicitly handled by the constructor here
}

#------------------------------------------------------------------------------#
# reIndexC.SparseBlockMatrixR
#  Re-indexing TO C for SparseBlockMatrixR objects
#
reIndexC.SparseBlockMatrixR <- function(x){
    #
    # Using lapply does NOT work here: if one of the list elements is an empty vector, adding 1 will
    #  mysteriously coerce it to a numeric vector (should be integer!). Not sure why this happens, but
    #  we shouldn't allow R to secretly change our data types without our permission. Use a for loop
    #  instead.
    #
    # UPDATE 05/13/14: Using '1L' keeps everything as an integer. This makes sense since '1' is a numeric
    #                   literal in R, while the corresponding literal for an integer is '1L'.
    #
    x$rows <- lapply(x$rows, function(x){ x - 1L})
    if(length(x$blocks) > 0) x$blocks <- lapply(x$blocks, function(x){ x - 1L})

    x$start <- 0

    x
} # END REINDEXC.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# reIndexR.SparseBlockMatrixR
#  Re-indexing TO R for SparseBlockMatrixR objects
#
reIndexR.SparseBlockMatrixR <- function(x){
    #
    # Using lapply does NOT work here: if one of the list elements is an empty vector, adding 1 will
    #  coerce it to a numeric vector (should be integer!). Using '1L' keeps everything as an integer.
    #  This makes sense since '1' is a numeric literal in R, while the corresponding literal for an
    #  integer is '1L'.
    #
    x$rows <- lapply(x$rows, function(x){ x + 1L})
    if(length(x$blocks) > 0) x$blocks <- lapply(x$blocks, function(x){ x + 1L})

    x$start <- 1

    x
} # END REINDEXR.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# SparseBlockMatrixR.list
#  List constructor
#
SparseBlockMatrixR.list <- function(x, ...){

    if( !is.list(x)){
        stop("Input must be a list!")
    }

    if( length(x) != 5 || !setequal(names(x), c("rows", "vals", "blocks", "sigmas", "start"))){
        stop("Input is not coercable to an object of type SparseBlockMatrixR, check list for the following elements: rows (list), vals (list), blocks (list), sigmas (numeric), start (integer)")
    }

    if(!is.list(x$rows) || !is.list(x$vals)){
        stop("rows and vals must both be lists of length pp!")
    }

    if( length(x$rows) != length(x$vals)){
        #
        # We enforce that rows and vals must have the same length, but relax this assumption for blocks and sigmas
        #  since they are mostly internal to the CCDr algorithm, and once the algorithm has run we may want to free
        #  up the memory associated with these elements
        #
        stop("rows and vals have different sizes; should all have the same length (pp)!!")
    }

    structure(x, class = "SparseBlockMatrixR")
} # END SPARSEBLOCKMATRIXR.LIST

#------------------------------------------------------------------------------#
# SparseBlockMatrixR.sparse
#  sparse object constructor
#
SparseBlockMatrixR.sparse <- function(x, sigmas, ...){

    if( !sparsebnUtils::is.sparse(x)){
        stop("Input must be a sparse object!")
    } else if(x$dim[1] != x$dim[2]){
        stop("Input must be square!")
    }

    pp <- x$dim[1]
    if(x$start == 0) x <- sparsebnUtils::reIndexR(x) # re-index rows and cols to start at 1 if necessary

    sbm.rows <- vector("list", length = pp)
    sbm.vals <- vector("list", length = pp)
    sbm.blocks <- vector("list", length = pp)

    if(missing(sigmas)){
        warning("Attempting to coerce sparse object to SparseBlockMatrixR with no data for sigmas: \n   Setting sigma_j = 0 by default.")
        sbm.sigmas <- rep(0, pp)
    } else{
        sbm.sigmas <- sigmas # just duplicating data here, but kept for consistency with other internal parameters for this method
    }


    # how to vectorize this???
    for(j in 1:pp){
        # Clear out the jth entry in the lists to be an empty vector
        sbm.rows[[j]] <- integer(0)
        sbm.vals[[j]] <- numeric(0)
        sbm.blocks[[j]] <- integer(0)
    }

    for(j in 1:pp){

        thisColIdx <- which(x$cols == j)
        rows <- as.integer(x$rows[thisColIdx])
        for(k in seq_along(rows)){
            row <- rows[k]

            sbm.rows[[j]] <- c(sbm.rows[[j]], row)
            sbm.rows[[row]] <- c(sbm.rows[[row]], j)

            sbm.vals[[j]] <- c(sbm.vals[[j]], x$vals[thisColIdx[k]])
            sbm.vals[[row]] <- c(sbm.vals[[row]], 0)

            # vals[rows[j][k]][block[j][k]] = beta_ji
            sbm.blocks[[j]] <- c(sbm.blocks[[j]], length(sbm.rows[[row]]))
            sbm.blocks[[row]] <- c(sbm.blocks[[row]], length(sbm.rows[[j]]))
        }

    }

    names(sbm.rows) <- names(sbm.vals) <- names(sbm.blocks) <- as.character(1:pp)

    #
    # NOTE: We use R-indexing by default. This can be changed by using reIndexC if necessary.
    #
    SparseBlockMatrixR(list(rows = sbm.rows, vals = sbm.vals, blocks = sbm.blocks, sigmas = sbm.sigmas, start = 1))
} # END SPARSEBLOCKMATRIXR.SPARSE

#------------------------------------------------------------------------------#
# SparseBlockMatrixR.matrix
#  matrix constructor
#
SparseBlockMatrixR.matrix <- function(x, sigmas, ...){

    if(nrow(x) != ncol(x)) stop("Input matrix must be square!")

    SparseBlockMatrixR(sparsebnUtils::as.sparse(x), sigmas, ...)
} # END SPARSEBLOCKMATRIXR.MATRIX

#------------------------------------------------------------------------------#
# as.list.SparseBlockMatrixR
#  Convert FROM SparseBlockMatrixR TO list
#  Even though internally the SBM object is a list, we must still manually define this function
#
as.list.SparseBlockMatrixR <- function(x){
    list(rows = x$rows, vals = x$vals, blocks = x$blocks, sigmas = x$sigmas, start = x$start)
} # END AS.LIST.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# as.matrix.SparseBlockMatrixR
#  Convert FROM SparseBlockMatrixR TO matrix
#
as.matrix.SparseBlockMatrixR <- function(x){
    pp <- length(x$rows)
    m <- matrix(0, nrow = pp, ncol = pp)

    ### 2015-03-02: Why was I using diag to construct this matrix?
    # m <- diag(rep(0, pp))

    if(x$start == 0) x <- sparsebnUtils::reIndexR(x)

    for(j in 1:pp){
        m[x$rows[[j]], j] <- x$vals[[j]]
    }

    ### 2015-03-02: Do not need to set dim attribute of matrix! (Already set by default constructor)
    # attributes(m)$dim <- c(pp, pp)
    # attributes(m)$dimnames <- list()
    rownames(m) <- as.character(1:nrow(m))
    colnames(m) <- as.character(1:ncol(m))

    m
} # END AS.MATRIX.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# edgeList.SparseBlockMatrixR
# Coerce SBM to edge list
#
#' @export
edgeList.SparseBlockMatrixR <- function(x){
    #
    # We have to be careful in obtaining the edge list of a SparseBlockMatrixR object:
    #  It is NOT the same as the rows slot since some of these components may have
    #  zero edge weights (see docs for SparseBlockMatrixR for explanation). Thus, in
    #  order to obtain the edge list, we need to check which indices in the rows slot
    #  have nonzero edge weights.
    #
    # y = rows, x = vals : Select the elements of rows which have nonzero values in vals,
    #                       accouting for possible round-off (hence .MACHINE_EPS).
    #
    el <- mapply(function(x, y){ y[which(abs(x) > sparsebnUtils::zero_threshold())]}, x$vals, x$rows)

    sparsebnUtils::edgeList(el)
} # EDGELIST.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# sparse.SparseBlockMatrixR
# 2016-01-22: Migrated to this file from s3-sparse.R
#
#' @export
sparse.SparseBlockMatrixR <- function(x, index = "R", ...){

    if(index != "R" && index != "C") stop("Invalid entry for index parameter: Must be either 'R' or 'C'!")

    pp <- length(x$rows)

    sp.rows <- integer(0)
    sp.cols <- integer(0)
    sp.vals <- numeric(0)

    sp.idx <- 0
    for(j in 1:pp){
        these.rows <- x$rows[[j]]
        these.vals <- x$vals[[j]]
        for(k in seq_along(these.rows)){

            # Only include nonzero values
            if(these.vals[k] != 0){
                sp.idx <- sp.idx + 1

                sp.rows <- c(sp.rows, these.rows[k])
                sp.cols <- c(sp.cols, j)
                sp.vals <- c(sp.vals, these.vals[k])
            }
        }
    }

    sp <- sparsebnUtils::sparse(list(rows = as.integer(sp.rows), cols = as.integer(sp.cols), vals = sp.vals, dim = c(pp, pp), start = 1))

    if(index == "R"){
        sp
    } else{
        sp$start <- 0
        sparsebnUtils::reIndexC(sp)
    }
} # END SPARSE.SPARSEBLOCKMATRIXR

# #------------------------------------------------------------------------------#
# # as.sparse.SparseBlockMatrixR
# #  Convert FROM SparseBlockMatrixR TO sparse
# #  By default, return the object using R indexing. If desired, the method can return C-style indexing by setting
# #    index = "C".
# # 2016-01-22: Migrated to this file from s3-sparse.R
# #
# as.sparse.SparseBlockMatrixR <- function(x, index = "R", ...){
#     sparse.SparseBlockMatrixR(x, index)
# } # END AS.SPARSE.SPARSEBLOCKMATRIXR

# to_graphNEL.SparseBlockMatrixR
#  Convert SBM object to graphNEL object
to_graphNEL.SparseBlockMatrixR <- function(x){
    ### This function require the 'graph' package to be installed
    if (!requireNamespace("graph", quietly = TRUE)) {
        stop("graph package (from BioConductor) required to coerce data to 'graphNEL' type!", call. = FALSE)
    }

    el <- sparsebnUtils::as.edgeList(x)
    el <- to_graphNEL(el)

    graph::graphNEL(nodes = as.character(1:num.nodes(x)), edgeL = el, edgemode = 'directed')
} # END TO_GRAPHNEL.SPARSEBLOCKMATRIXR

get.adjacency.matrix.SparseBlockMatrixR <- function(x){
    sparsebnUtils::get.adjacency.matrix(sparsebnUtils::as.edgeList(x))
} # END GET.ADJACENCY.MATRIX.SPARSEBLOCKMATRIXR

num.nodes.SparseBlockMatrixR <- function(x){
    ### The number of nodes should be exactly the same as the length of the rows list
    length(x$rows)
} # END NUM.NODES.SPARSEBLOCKMATRIXR

num.edges.SparseBlockMatrixR <- function(x){
    ### The number of nodes should be exactly the same as the length of the rows list
    sparsebnUtils::num.edges(sparsebnUtils::as.edgeList(x))
} # END NUM.EDGES.SPARSEBLOCKMATRIXR

# This function is (so far) only used in unit tests
is.zero.SparseBlockMatrixR <- function(x){
    check_if_zero <- (length(unlist(x$x$rows)) == 0)

    check_if_zero
} # END IS.ZERO.SPARSEBLOCKMATRIXR

#------------------------------------------------------------------------------#
# .init_sbm
# Internal function for initializing a SparseBlockMatrixR object directly
#  from a matrix AND a sigmas vector
#
.init_sbm <- function(init_matrix, init_sigmas){
    stopifnot(sparsebnUtils::check_if_matrix(init_matrix))
    stopifnot(nrow(init_matrix) == ncol(init_matrix))

    stopifnot(is.numeric(init_sigmas))
    stopifnot(length(init_sigmas) == nrow(init_matrix))

    sbm <- SparseBlockMatrixR(init_matrix, init_sigmas)
    ### 2016-04-29 These two lines replaced with the above call
    # sbm <- suppressWarnings(SparseBlockMatrixR(init_matrix)) # suppress warnings since we are working in a controlled environment
    # sbm$sigmas <- init_sigmas

    sbm
} # END .INIT_SBM

#------------------------------------------------------------------------------#
# .to_B.SparseBlockMatrixR
# Internal function to convert estimates from the (Rho, R) parametrization to
#  the standard (B, Omega) parametrization.
#
to_B.SparseBlockMatrixR <- function(sbm){
    ### Need to re-parametrize the betas FIRST
    sbm$vals <- mapply(function(x, y) x/y, sbm$vals, sbm$sigmas) # Divide each vals vector by sigmas[j]
    sbm$sigmas <- 1/ (sbm$sigmas)^2

    sbm
}
