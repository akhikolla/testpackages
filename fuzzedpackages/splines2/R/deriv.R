##
## R package splines2 by Wenjie Wang and Jun Yan
## Copyright (C) 2016-2020
##
## This file is part of the R package splines2.
##
## The R package splines2 is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package splines2 is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##

##' Derivatives of Spline Bases
##'
##' Returns derivatives of given order for the spline bases.
##'
##' At knots, the derivative is defined to be the right derivative. By default,
##' the function returns the first derivatives.  For derivatives of order
##' greater than one, the nested call such as \code{deriv(deriv(expr))} is
##' supported but not recommended.  For a better performance, argument
##' \code{derivs} should be specified instead.
##'
##' This function is designed for objects produced by this package.  It
##' internally extracts necessary specification about the spline/polynomial
##' basis matrix from its attributes. Therefore, the function will not work if
##' the key attributions are not available after some operations.
##'
##' @name deriv
##'
##' @param expr Objects of class \code{bSpline2}, \code{ibs}, \code{dbs},
##'     \code{mSpline}, \code{iSpline}, or \code{cSpline}, etc.
##' @param derivs A positive integer specifying the order of derivatives. By
##'     default, it is \code{1L} for the first derivative.
##' @param ... Other arguments that are not used.
##'
##' @return A numeric matrix of the same dimension with the input \code{expr}.
##'
##' @example inst/examples/ex-deriv.R
##'
##' @importFrom stats deriv
NULL


##' @rdname deriv
##' @export
deriv.bSpline2 <- function(expr, derivs = 1L, ...)
{
    attr(expr, "derivs") <- derivs
    ## checks if key attributions still exist
    check_attr(expr, check_derivs = FALSE)
    dMat <- do.call(dbs, attributes(expr))
    class(dMat) <- c("matrix", "dbs")
    dMat
}


##' @rdname deriv
##' @export
deriv.dbs <- function(expr, derivs = 1L, ...)
{
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    ## checks if key attributions still exist
    check_attr(expr, check_derivs = TRUE)
    dMat <- do.call(dbs, attributes(expr))
    class(dMat) <- c("matrix", "dbs")
    dMat
}


##' @rdname deriv
##' @export
deriv.ibs <- function(expr, derivs = 1L, ...)
{
    ## quick check on derivs
    derivs <- as.integer(derivs)
    if (derivs < 1L) {
        stop("'derivs' has to be a positive integer.")
    }
    ## checks if key attributions still exist
    check_attr(expr, check_derivs = FALSE)
    ## if first derivative, take result from existing attribute if exists
    if (derivs == 1L) {
        out <- attr(expr, "bsMat")
        if (is.null(out))
            out <- do.call(bSpline, attributes(expr))
        return(out)
    }
    ## for higher order of derivative
    attr(expr, "derivs") <- derivs - 1L
    dMat <- do.call(dbs, attributes(expr))
    class(dMat) <- c("matrix", "dbs")
    dMat
}


##' @rdname deriv
##' @export
deriv.mSpline <- function(expr, derivs = 1L, ...)
{
    derivs0 <- attr(expr, "derivs")
    ## checks if key attributions still exist
    check_attr(expr, check_derivs = FALSE)
    attr(expr, "derivs") <- ifelse(is.null(derivs0), derivs, derivs0 + derivs)
    dMat <- do.call(mSpline, attributes(expr))
    class(dMat) <- c("matrix", "mSpline")
    dMat
}


##' @rdname deriv
##' @export
deriv.iSpline <- function(expr, derivs = 1L, ...)
{
    ## checks if key attributions still exist
    check_attr(expr, check_derivs = FALSE)
    attr(expr, "derivs") <- derivs - 1L
    dMat <- do.call(mSpline, attributes(expr))
    class(dMat) <- c("matrix", "mSpline")
    dMat
}


##' @rdname deriv
##' @export
deriv.cSpline <- function(expr, derivs = 1L, ...)
{
    ## checks if key attributions still exist
    check_attr(expr, check_derivs = FALSE)
    scl <- attr(expr, "scales")
    ## not scaled
    if (is.null(scl)) {
        derivs <- as.integer(derivs)
        if (derivs == 1L) {
            return(do.call(iSpline, attributes(expr)))
        }
        if (derivs == 2L) {
            return(do.call(mSpline, attributes(expr)))
        }
        attr(expr, "derivs") <- derivs - 2L
        return(do.call(mSpline, attributes(expr)))
    }
    ## else scaled
    derivs0 <- attr(expr, "derivs")
    attr(expr, "derivs") <- ifelse(is.null(derivs0), derivs, derivs0 + derivs)
    do.call(cSpline, attributes(expr))
}


##' @rdname deriv
##' @export
deriv.bernsteinPoly <- function(expr, derivs = 1L, ...)
{
    ## checks if key attributions still exist
    check_attr(expr,
               attrs = c("x", "degree", "Boundary.knots",
                         "intercept", "integral"))
    attr(expr, "derivs") <- attr(expr, "derivs") + derivs
    do.call(bernsteinPoly, attributes(expr))
}
