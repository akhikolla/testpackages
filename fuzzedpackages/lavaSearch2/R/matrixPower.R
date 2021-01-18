### matrixPower.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 23 2017 (16:52) 
## Version: 
## last-updated: nov  2 2018 (15:00) 
##           By: Brice Ozenne
##     Update #: 59
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * matrixPower (documentation)
##' @title Power of a Matrix
##' @description Compute the power of a matrix.
##' @name matrixPower
##' 
##' @param object a matrix.
##' @param power [numeric] power to be applied to the matrix.
##' @param symmetric [logical] is the matrix symmetric? Argument passed to the function \code{eigen}.
##' @param tol [numeric >0] the threshold under which the eigenvalues are set to 0.
##' @param print.warning [logical] should a warning be print when some or the eigenvalues are not strictly positive.
##' 
##' @return A matrix.

## * matrixPower (example)
##' @rdname matrixPower
##' @examples
##' ## symmetric matrix
##' set.seed(10)
##' M <- matrix(rnorm(20*6),20,6)
##' Sigma <- var(M)
##' Sigma.half <- matrixPower(Sigma, power = 1/2, symmetric = TRUE)
##' round(Sigma.half %*% Sigma.half - Sigma,5)
##' 
##' iSigma <- matrixPower(Sigma, power = -1, symmetric = TRUE)
##' round(iSigma %*% Sigma,5)
##' 
##' iSigma.half <- matrixPower(Sigma, power = -1/2, symmetric = TRUE)
##' round(iSigma.half %*% iSigma.half - iSigma,5)
##' 
##' ## non symmetric matrix
##' set.seed(10)
##' M <- matrix(abs(rnorm(9)), 3, 3) + diag(1,3,3)
##' M-t(M)
##' 
##' iM <- matrixPower(M, power = -1, symmetric = FALSE)
##' round(iM %*% M,5)
##' 
##' iM.half <- matrixPower(M, power = -1/2, symmetric = FALSE)
##' round(iM.half %*% iM.half %*% M,5)
##' 

## * matrixPower (code)
##' @rdname matrixPower
##' @export
matrixPower <- function(object, power, symmetric, tol = 1e-12, print.warning = TRUE){
    object.eigen <- eigen(object, symmetric = symmetric)

    warn <- 0
    if(power<0 && any(object.eigen$values<tol)){
        if(print.warning){
            warning("small or negative eigenvalues are set to ",tol,"\n")
        }
        warn <- 1
        object.eigen$values[object.eigen$values < tol] <- tol
    }
    nRow <- length(object.eigen$values)
    D <- diag(object.eigen$values^power, nrow = nRow, ncol = nRow)

    if(symmetric){
        out <- object.eigen$vectors %*% D %*% t(object.eigen$vectors)
    }else{
        out <- object.eigen$vectors %*% D %*% solve(object.eigen$vectors)
    }

    if(warn==1){
        attr(out,"warning") <- "small or negative eigenvalues"
    }
    return(out)

}

#----------------------------------------------------------------------
### matrixPower.R ends here
