
## DO NOT ROXYGENISE THAT!!!
##
##' Check the gradient provided in a \code{covMan} object.
##'
##' Each of the two matrices \code{x1} and \code{x2} with \code{n1}
##' and \code{n2} rows can be given or instead be drawn at random. The
##' matrix of kernel values with dimension \code{c(n1, n2)} is
##' computed, together with its gradient with dimension \code{c(n1,
##' n2, npar)} where \code{npar} is the number of parameters of the
##' kernel. A numerical differentiation w.r.t. the kernel parameters
##' is performed for the kernel value at \code{x1} and \code{x2}, and
##' the result is compared to that provided by the kernel function
##' (the function described in the slot named \code{"kernel"} of
##' \code{object}). Note that the value of the parameter vector is the
##' value provided by \code{coef(object)} and it can be changed by
##' using the replacement method \code{`coef<-`} if needed.
##' 
##' @title Check the gradient provided in a \code{covMan} object
##'
##' @param object A \code{covMan} object.
##'
##' @param sym Logical. If \code{TRUE}, the check is done assuming
##' that \code{x2} is identical to \code{x1}, so the provided values
##' for \code{x2} and \code{n2} (if any) will be ignored.
##'
##' @param x1 Matrix to be used as the first argument of the kernel.
##' 
##' @param n1 Number of rows for the matrix \code{x1}. Used only when
##' \code{x1} is not provided.
##'
##' @param x2 Matrix to be used as the second argument of the kernel.
##' 
##' @param n2 Number of rows for the matrix \code{x2}. Used only when
##' \code{x2} is not provided.
##'
##' @param XLower Vector of lower bounds to draw \code{x1} and
##' \code{x2} when needed.
##'
##' @param XUpper Vector of upper bounds to draw \code{x1} and
##' \code{x2} when needed.
##'
##' @param plot \code{Logical}. If \code{TRUE}, a plot is shown
##' comparing the two arrays of gradients.
##' 
##' @return A list of results related to the Jacobians
##'
##' \itemize{
##'
##' \item{\code{test}}{
##' Max of the absolute difference between the gradient obtained by
##' numeric differentiation and the gradient provided by the kernel
##' object.
##' }
##' 
##' \item{\code{Jnum}, \code{J}}{
##' Jacobians (arrays) computed with \code{numDeriv::jacobian} and
##' provided by the kernel object.
##' }
##'
##' \item{\code{x1}, \code{x2}, \code{K}}{ The matrices used for the
##' check, and the matrix of kernel values with dimension \code{c(n1,
##' n2)}. The element \code{x2} can be \code{NULL} if the
##' determination of the matrix \code{x2} was not necessary. }
##'
##' }
##'
##' @note
##'
##' As a rule of thumb, a gradient coded without error gives a value
##' of \code{test} less than \code{1e-4}, and usually the value is
##' much smaller than that.
##'
##' @section Caution: For now the function only works when
##' \code{object} has class \code{"covMan"}.
##' 
##' @author Yves Deville
##' 
checkGrad <- function(object,
                      sym = TRUE,
                      x1 = NULL, n1 = 10, 
                      x2 = NULL, n2 = NULL, 
                      XLower = NULL, XUpper = NULL,
                      plot = TRUE) {

    if (class(object) != "covMan") {
        stop("For now this function works only when 'object'\n",
             "has class \"covMan\"") 
    }

    if (!object@acceptMatrix) {
        if (!sym) {
            warning("when 'object' has its slot \"acceptMatrix\"\n",
                    "FALSE, 'sym' can only be TRUE. Forced.")
        }
    }
    
    if (!object@hasGrad) 
        stop("'object' does not compute the gradient, so no check is possible")
    
    d <- object@d
    
    ##=========================================================================
    ## Prepare matrices 'x1' and 'x2' if needed.
    ##=========================================================================
    
    if (is.null(x1) || (!sym && is.null(x2))) {

        if (is.null(XLower)) {
            XLower <- rep(0.0, d)
        } else if (length(XLower) != d) {
            stop("'XLower' must have length ", d)
        }
        
        if (is.null(XUpper)) {
            XUpper <- rep(1.0, d)
        } else {
            if (length(XUpper) != d) {
                stop("'XUpper' must have length ", d)
            }
            if (any(XUpper <= XLower)) {
                stop("'XUpper' must be elementwise greater than 'Xlower'")
            }
        }
        XRange <- XUpper - XLower
    }
    
    if (is.null(x1)) {
        x1 <- array(runif(n1 * d), dim = c(n1, d))
        x1 <- sweep(x1, MARGIN = 2, STATS = XRange, FUN = "*")
        x1 <- sweep(x1, MARGIN = 2, STATS = XLower, FUN = "+")
        colnames(x1) <- inputNames(object)
    } 
    
    if (sym) {
        if (!is.null(x2) || !is.null(n2)) {
            warning("since 'sym' is TRUE, the provided value of 'x2' or 'n2' is ignored")
        }
        ## if the kernel function does not accept 'x2 = NULL'
        if (!is.null(formals(object@kernel)$x2)) {
            x2 <- x1
        }
        n2 <- n1
    } else if (is.null(x2)) {
        x2 <- array(runif(n2 * d), dim = c(n2, d))
        x2 <- sweep(x2, MARGIN = 2, STATS = XRange, FUN = "*")
        x2 <- sweep(x2, MARGIN = 2, STATS = XLower, FUN = "+")
        colnames(x2) <- inputNames(object)
    }
    
    ##=========================================================================
    ## Compute the jacobian with the 'numDeriv' package.
    ## 
    ## This could be achieved by copying the kernel fun and then
    ## hacking the formals as in
    ##
    ## formals(KFun) <- alist(par = , x1 = , x2 = )
    ##
    ##=========================================================================
    
    cpar <- coef(object)

    if (object@acceptMatrix) {

        KFun <- function(par, x1, x2) {
            object@kernel(x1 = x1, x2 = x2, par = par)
        }
        
        K1 <- KFun(cpar, x1 = x1, x2 = x2)   
        Jnum <- jacobian(func = KFun, x = cpar, x1 = x1, x2 = x2)
    
        ## reshape as an array with suitable dims.
        J2 <- array(Jnum, dim = c(n1, n2, object@parN),
                    dimnames = list(rownames(x1), rownames(x2), object@kernParNames))
        for (j in 1L:object@parN) J2[ , , j] <- Jnum[ , j]

        
        J <- attr(K1, "gradient")
        attr(K1, "gradient") <- NULL
    
        ## If 'J' is a list (hence not an array) coerce it manually
        if (is.list(J)) {
            J1 <- array(NA, dim = c(n1, n2, object@parN),
                        dimnames = list(rownames(x1), rownames(x2), object@kernParNames))
            for (j in 1L:object@parN) J1[ , , j] <- J[[j]]
            J <- J1
        }
        
        if (plot) {
            plot(J2, pch = 21, col = "SpringGreen3", cex = 0.6, ylim = range(J, J2), 
                 main = c("orange: kernel, green: num. diff."))
            points(J, pch = 16, col = "orangered", cex = 0.4)
        }
        
    }  else {
        
        KFunNoMat <- function(par, x1, compGrad, index = 0) {
            newobj <- object
            coef(newobj) <- par
            covMat(newobj, X = x1, compGrad = compGrad, index = index,
                   checkNames = FALSE)
        }
        
        K1 <- KFunNoMat(cpar, x1 = x1, compGrad = FALSE)   
        Jnum <- jacobian(func = KFunNoMat, x = cpar, x1 = x1, compGrad = FALSE)
        
        ## reshape as an array with suitable dims.
        J2 <- array(Jnum, dim = c(n1, n1, object@parN),
                    dimnames = list(rownames(x1), rownames(x1), object@kernParNames))
        for (j in 1L:object@parN) J2[ , , j] <- Jnum[ , j]
        
        J <- array(NA, dim = c(n1, n1, object@parN),
                   dimnames = list(rownames(x1), rownames(x1), object@kernParNames))
        
        ## reshape as an array with suitable dims.
        for (j in 1:object@parN) {
            Kprov <- KFunNoMat(cpar, x1 = x1, compGrad = TRUE, index = j)   
            J[ , , j] <- attr(Kprov, "gradient")
        }

        
        if (plot) {
            plot(J2, pch = 21, col = "SpringGreen3", cex = 0.6, ylim = range(J, J2), 
                 main = c("orange: kernel, green: num. diff."))
            points(J, pch = 16, col = "orangered", cex = 0.4)
        }
    }
    
    list(test = max(abs(J2 - J)),
         Jnum = J2, J = J,
         x1 = x1, x2 = x2, K = K1)
              
}
