#
#  ccdrAlgorithm-functions.R
#  ccdrAlgorithm
#
#  Created by Dacheng Zhang on 9/3/16.
#  Copyright (c) 2014-2017 Bryon Aragam. All rights reserved.
#

## returns TRUE if ivn_list is a list of vectors or NULL elements,
check_if_ivn_list <- function(ivn) {
    ## check if it is a list
    if(!is.list(ivn)) return(FALSE)

    ## check if every component is a vector of NULL
    return(all(sapply(ivn, is.vector) | sapply(ivn, is.null)))
} # END CHECK_IF_IVN_LIST

## returns TRUE if ivn_list has length nn, the number of sample rows
check_ivn_size <- function(ivn, data) {
    ## check if length matches with nn
    return(length(ivn) == nrow(data))
} # END CHECK_IF_IVN_SIZE

## returns TRUE if a vector component of 'ivn' is NULL,
## or has all correct labels of nodes under intervention in this sample:
## 1) integer, 2) between 1 and pp, and 3) no duplicates
check_vector_label <- function(vec, pp) {

    if(is.null(vec)) return(TRUE)

    ## Note: If a vector has only integers and NAs, is.integer returns all TRUE
    ## e.g.: c(NA, 1L, NA, 3L, NA, 5L)
    ## However, c(1L, NA, 3L, 4, NA) returns all FALSE
    ## check if labels are integers
    if(any(is.na(vec)) || !is.numeric(vec)) {
        stop("Non-integer label(s) found in one or more components in ivn.")
        return(FALSE)
    }

    ## check if labels are in 1..pp
    if(any(vec < 1) | any(vec > pp)) {
        stop(sprintf("Labels should all be between 1 and %d to refer to the columns of data.", pp))
        return(FALSE)
    }

    ## check if labels are unique
    if(anyDuplicated(vec)) {
        stop("Duplicated label(s) found in one component in ivn.")
        return(FALSE)
    }

    return(TRUE)
} # END CHECK_VECTOR_LABEL

## returns TRUE if every vector in 'ivn' is NULL,
## or has correct labels: integer, between 1 and pp, and no duplicates
check_ivn_label <- function(ivn, data) {
    sapply(ivn, check_vector_label, ncol(data))
} # END CHECK_IVN_LABEL

