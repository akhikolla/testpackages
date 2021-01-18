### symmetrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:42) 
## Version: 
## Last-Updated: feb 19 2018 (09:33) 
##           By: Brice Ozenne
##     Update #: 19
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Symmetrize a Matrix
#' @description Complete the upper (or lower) extra-diagonal terms  in order to obtain a symmetric matrix.
#' @param M a matrix.
#' @param update.upper [logical] should the upper extra diagonal terms be updated using the lower extra diagonal terms?
#' 
#' @examples
#' symmetrize <- lavaSearch2:::symmetrize
#'
#' ## example
#' M <- matrix(NA, 4, 4)
#' M[lower.tri(M)] <- 1:6
#'
#' symmetrize(M, update.upper = TRUE) # good
#'
#' M[upper.tri(M, diag = FALSE)] <- M[lower.tri(M, diag = FALSE)]
#' M # wrong
#' @keywords internal
symmetrize <- function(M, update.upper = TRUE){

    if(is.null(update.upper)){
        M <- (M+t(M))
        diag(M) <- diag(M)/2
    }else if(update.upper){
        M[upper.tri(M, diag = FALSE)] <- t(M)[upper.tri(M, diag = FALSE)]
    }else{
        M[lower.tri(M, diag = FALSE)] <- t(M)[lower.tri(M, diag = FALSE)]
    }
    return(M)    
}

##----------------------------------------------------------------------
### symmetrize.R ends here
