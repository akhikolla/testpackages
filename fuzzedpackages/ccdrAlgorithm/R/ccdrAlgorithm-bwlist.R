#
#  ccdrAlgorithm-bwlist.R
#  ccdrAlgorithm
#
#  Created by Bryon Aragam (local) on 8/11/17.
#  Copyright (c) 2014-2017 Bryon Aragam. All rights reserved.
#

#
# PACKAGE CCDRALGORITHM: Helper methods for black/white lists
#
#   CONTENTS:
#     names_to_indices
#     rows_to_list
#     bwlist_check
#     bwlist_to_weights
#

### Just a wrapper for match with a better name
names_to_indices <- function(v, names){
    match(v, names)
} # END NAMES_TO_INDICES

### Returns a list whose components are the rows of a matrix
rows_to_list <- function(m){
    lapply(1:nrow(m), function(j) m[j,])
} # END ROWS_TO_LIST

### Check correctness of input and transform from matrix to list of indices
bwlist_check <- function(bwlist, names){
    ## Consistency checks
    if(!is.matrix(bwlist) || ncol(bwlist) != 2){
        stop("Input must be a matrix with exactly 2 columns!")
    }

    if(any(is.na(bwlist))){
        stop("Input cannot have missing values!")
    }

    ### Convert characters names to indices
    if(is.character(bwlist)){
        bwlist <- as.vector(bwlist)
        bwlist <- names_to_indices(bwlist, names)
        bwlist <- matrix(bwlist, ncol = 2)
    }

    storage.mode(bwlist) <- "integer" # This is important in ccdr_call to check overlap between blacklist and whitelist, fails if numerics are mixed with ints
    rows_to_list(bwlist)
} # END BWLIST_CHECK

### Convert b/w lists to weight matrix of {-1,0,1}
#     -1 = black listed (guaranteed to be absent / zero)
#      0 = white listed (guaranteed to be present / nonzero)
#      1 = gray listed (may or may not be final model)
bwlist_to_weights <- function(black, white, nnode){
    weights <- matrix(1L, ncol = nnode, nrow = nnode)

    if(!is.null(white)){
        for(k in 1:length(white)){
            weights[white[[k]][1], white[[k]][2]] <- 0L
        }
    }

    if(!is.null(black)){
        for(k in 1:length(black)){
            weights[black[[k]][1], black[[k]][2]] <- -1L
        }
    }

    weights
} # END BWLIST_TO_WEIGHTS
