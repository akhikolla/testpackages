##  DO NOT ROXYGENIZE THIS! This is a skeleton for the .Rd
##  file that might have been edited later because 'kergp' does not
##  Roxygenize its man.
##=============================================================================
##' Simulation of paths from a \code{gp} object.
##'
##' @title Simulation of paths from a \code{gp} object.
##' 
##' @param object An object with class \code{"gp"}.
##' 
##' @param nsim Number of paths wanted.
##'
##' @param seed Not used yet.
##'
##' @param newdata A data frame containing the inputs values used for
##' simulation as well as the required trend covariates, if any. This
##' is similar to the \code{newdata} formal in
##' \code{\link{predict.gp}}.
##'
##' @param cond Logical. Should the simulations be conditional on
##' the observations used in the object or not?
##' 
##' @param trendKnown Logical. If \code{TRUE} the vector of trend
##' coefficients will be regarded as known so all simulated paths
##' share the same trend. When \code{TRUE}, the trend must have been
##' estimated so that its estimation covariance is known. Then each
##' path will have a trend
##'
##' @param newVarNoise Variance of the noise for the "new" simulated
##' observations. For the default \code{NULL}, the noise variance
##' found in \code{object} is used. Note that if a very small positive
##' value is used, each simulated path is the sum of the trend the
##' smooth GP part and an interval containing say \eqn{95}\% of the
##' simultated responses can be regarded as a confidence interval
##' rather than a prediction interval.
##' 
##' @param output The type of output wanted. A simple matrix as in
##' standard simulation methods may be quite poor, since interesting
##' intermediate results are then lost.
##'
##' @param label A label that will be attached to the object if
##' \code{output} is \code{"list"}.
##'
##' @param unit A unit that will be attached to the object if
##' \code{output} is \code{"list"}.
##' 
##' @return A matrix with the simulated paths as its columns or a more
##' complete list with more results. This list which is given the S3 class
##' \code{"simulate.gp"} has the following elements.
##'
##' \itemize{
##' \item{\code{X}, \code{F}, \code{y}}{
##' Inputs, trend covariates and response.
##' }
##' \item{\code{XNew}, \code{FNew}}{
##' New inputs, new trend covariates.
##' }
##' \item{\code{sim}}{
##' Matrix of simulated paths.
##' }
##' \item{\code{trend}}{
##' Matrix of simulated trends.
##' }
##' \item{\code{trendKnown}, \code{noise}, \code{newVarNoise}}{
##' Values of the formals.
##' }
##' \item{\code{Call}}{
##' The call.
##' }
##' }
##'
##' 
##' @note When \code{betaKnown} is \code{FALSE}, the simulated
##' \emph{trend} and the \emph{smooth GP} parts of a simulation are
##' usually correlated, and their sum will show less dipersion than
##' each of the two components. The covariance of the vector
##' \eqn{\widehat{\boldsymbol{\beta}}}{\beta} can be regarded as the
##' posterior distribution corresponding to a non-informative prior,
##' and the distribution from which a new path is drawn is the
##' predictive distribution.
##' 
##' @author Yves Deville
##'
##' @examples
##' set.seed(314159)
##' n <- 40
##' x <- sort(runif(n))
##' y <- 2 + 4 * x  + 2 * x^2 + 3 * sin(6 * pi * x ) + 1.0 * rnorm(n)
##' nNew <- 60; xNew <- sort(runif(nNew))
##' df <- data.frame(x = x, y = y)
##'
##' ##=========================================================================
##' ## use a Matern 3/2 covariance
##' ##=========================================================================
##' myKern <- k1Matern3_2
##' inputNames(myKern) <- "x"
##' mygp1 <- gp(formula = y ~ x + I(x^2) + sin(6 * pi * x),
##'             data = df, inputs = "x",
##'             parCovLower = c(0.01, 0.01), parCovUpper = c(10, 100),
##'             cov = myKern, estim = TRUE, noise = TRUE)
##' mygp2 <- gp(formula = y ~ sin(6 * pi * x),
##'             data = df, inputs = "x",
##'             parCovLower = c(0.01, 0.01), parCovUpper = c(10, 100),
##'             cov = myKern, estim = TRUE, noise = TRUE)
##'
##' ##=========================================================================
##' ## New data
##' ##=========================================================================
##' nNew <- 300
##' xNew <- seq(from = -0.2, to= 1.2, length.out = nNew)
##' dfNew <- data.frame(x = xNew)
##' 
simulate.gp <- function(object, nsim = 1L, seed = NULL,
                        newdata = NULL,
                        cond = TRUE,
                        trendKnown = FALSE,
                        newVarNoise = NULL,
                        checkNames = TRUE,
                        output = c("list", "matrix"),
                        label = "y", unit = "",
                        ...) {
    
    mc <- match.call()
    output <- match.arg(output)
    noise <- object$noise
    
    if (is.null(newVarNoise)) {
        if (noise) {
            newVarNoise <- object$varNoise
        } 
    } else {
        if (newVarNoise < 0.0) stop("'varNoise' must be >= 0.0")
        noise <- (newVarNoise > 0.0) 
    }
   
    N <- nsim
    
    n <- nrow(object$F)
    p <- ncol(object$F)
    d <- ncol(object$X)
    nms <- object$inputNames
    
    if (missing(newdata)) {
        XNew <- object$X
        FNew <- object$F
    } else {
        XNew <- newdata[ , nms, drop = FALSE]
        nms <- object$inputNames   
        New <- newdata[ , nms, drop = FALSE]
        tt <- delete.response(terms(object))
        mf <- model.frame(tt, data = data.frame(newdata))
        FNew <- model.matrix(tt, data = mf)
    }
    nNew <- nrow(XNew)
    
    if (!cond) {

        if (!trendKnown) {
            ## Simulate 'beta'...
            Znorm <- array(rnorm(p * N), dim = c(p, N))
            dBetaSim <- backsolve(object$RStar, Znorm)
            trendNewSim <-  FNew %*% sweep(x = dBetaSim, MARGIN = 1,
                                           STATS = object$betaHat, FUN = "+")
        }  else {
            ## ... or not.
            trend <- FNew %*% object$betaHat
            trendNewSim <-  array(trend, dim = c(nNew, N))
        }

        ## simulate the smooth GP part
        zetaNewSim <- simulate(object = object$covariance,
                               nsim = nsim,
                               seed = seed,
                               X = XNew, checkNames = checkNames, ...)
        
        yNewSim <- trendNewSim + zetaNewSim

        ## add noise if needed
        if (noise) {
            yNewSim <- yNewSim + sqrt(newVarNoise) *
                array(rnorm(nNew * N), dim = c(nNew, N)) 
        }
        
        if (output == "list") {
            yNewSim <- list(X = object$X, F = object$F, y = object$y,
                            XNew = XNew, FNew = FNew,
                            sim = yNewSim,
                            trend = trendNewSim,
                            trendKnown = trendKnown,
                            noise = noise,
                            newVarNoise = newVarNoise,
                            label = label,
                            unit = unit,
                            Call = mc)
            class(yNewSim) <- "simulate.gp"
            return(yNewSim)   
        } else {
            attr(yNewSim, "trendKnown") <- trendKnown
            return(yNewSim)
        }
   
    }
    
    ##===========================================================================
    ## compute 'KCond' and its Cholesky root as needed in step 2. 'K'
    ## could be retreived from its Cholesky root?
    ##=========================================================================== 
    K <- covMat(object$covariance, X = object$X, checkNames = checkNames,
                compGrad = FALSE)
    KStar <- forwardsolve(object$L, K)
    KCond <- K - crossprod(KStar)
    LCond <- t(chol(KCond))

    if (!trendKnown) {
    
        ##=======================================================================
        ## STEP 1 simulate 'beta - betaHat'.  'dBetaSim' is a matrix
        ## with random (betaSim - betaHat) as its columns
        ## ======================================================================
        Znorm <- array(rnorm(p * N), dim = c(p, N))
        dBetaSim <- backsolve(object$RStar, Znorm)
        trendNewSim <-  FNew %*% sweep(x = dBetaSim, MARGIN = 1,
                                       STATS = object$betaHat, FUN = "+")
        
        ##=======================================================================
        ## STEP 2 simulate 'zeta' conditional on 'beta' and 'y'
        ##=======================================================================
    
        ## 'dBetaSim' -> the 2nd part of L^{-1} [y - F betaSim]
        E <- -object$FStar %*% dBetaSim
        
        ## add the 1st part of L^{-1} [y - F betaSim], i.e.  'eStar'
        E <- sweep(x = E, MARGIN = 1, STATS = object$eStar, FUN = "+")
    
        ## left mutliply by t(KStar) = K %*% L^{-T}
        E <- t(KStar) %*% E
        
        ## now, find the 'zetaSim' conditional on 'y' and 'betaSim' 
        ZSim <- array(rnorm(n * N), dim = c(n, N))
        ZSim <- E + LCond %*% ZSim
        
    } else {
        trend <- FNew %*% object$betaHat
        trendNewSim <-  array(trend, dim = c(nNew, N))
        ## now, find the 'zetaSim' conditional of 'y' with 'beta' fixed
        ZSim <- array(rnorm(n * N), dim = c(n, N))
        ZSim <- LCond %*% ZSim
        eStarMod <- t(KStar) %*% object$eStar
        ZSim <- sweep(x = ZSim, MARGIN = 1, STATS = eStarMod, FUN = "+")
    }
    
    ##=======================================================================
    ## STEP 3: now krige the N columns to find 'ZNewSim'
    ##=======================================================================
    
    kNew <- covMat(object$covariance, X = object$X, Xnew = XNew, compGrad = FALSE)
    KNew <- covMat(object$covariance, X = XNew, Xnew = XNew, compGrad = FALSE)
    kNewStar <- forwardsolve(object$L, kNew)
    
    KNewCond <- KNew - crossprod(kNewStar)
    LNewCond <- t(chol(KNewCond))
    
    ZNewSim <- array(rnorm(nNew * N), dim = c(nNew, N))
    ZNewSim <- LNewCond %*% ZNewSim
    ZNewSim <- ZNewSim + t(kNewStar) %*% forwardsolve(object$L, ZSim)

    yNewSim <- trendNewSim + ZNewSim
    
    if (noise) {
        yNewSim <- yNewSim + sqrt(newVarNoise) *
            array(rnorm(nNew * N), dim = c(nNew, N)) 
    }
        
    if (output == "list") {
        yNewSim <- list(X = object$X, F = object$F, y = object$y,
                        XNew = XNew, FNew = FNew,
                        sim = yNewSim,
                        trend = trendNewSim,
                        trendKnown = trendKnown,
                        noise = noise,
                        newVarNoise = newVarNoise,
                        label = label,
                        unit = unit,
                        Call = mc)
        class(yNewSim) <- "simulate.gp"
        return(yNewSim)   
    } else {
        attr(yNewSim, "trendKnown") <- trendKnown
        return(yNewSim)
    }
        
}
##' Function to plot simulations form a \code{gp} object.
##'
##' @title  Function to plot simulations form a \code{gp} object
##'
##' @param x An object containing simulations, produced by 'simulate'
##' with \code{output = "list"}.
##'
##' @param y Not used yet.
##'
##' @param col Named list of colors to be used, with elements
##' \code{"sim"} and \code{"trend"}.
##' 
##' @param show A logical vector telling which elements must be shown.
##'
##' @param ... Further argument passed to \code{plot}.
##' 
##' @return Nothing
##'
##' @note For now, this function can be used only when the number of
##' inputs is one.
##'
##' @seealso \code{\link{simulate.gp}}.
plot.simulate.gp <- function(x, y,
                             col = list("sim" = "SpringGreen3", "trend" = "orangered"),
                             show = c("sim" = TRUE, "trend" = TRUE, "y" = TRUE),
                             ...) {
    
    if (ncol(x$X) > 1L) {
        stop("For now, this function only for the one-dimensional case")
    }
    
    X <- drop(as.matrix(x$X))
    XNew <- drop(as.matrix(x$XNew))

    ylab <- x$label
    if (nchar(x$unit)) ylab <- sprintf("%s (%s)", ylab, x$unit)
    
    plot(X, x$y, type = "n",
         xlim = range(X, XNew),
         ylim = range(x$sim, x$trend, x$y),
         ylab = ylab,
         ...)


    ## find a value of the opacitiy level 'alpha'
    nSim <- ncol(x$sim)
    alpha <- exp(-sqrt(((nSim - 1)/ 200)))
    
    if (show["trend"]) {
        matlines(XNew, x$trend, pch = 16,
                 col = translude(col[["trend"]], alpha = alpha))
    }
    if (show["sim"]) {
        matlines(XNew, x$sim, pch = 16,
                 col = translude(col[["sim"]], alpha = alpha))
    }
    if (show["y"]) {
        points(X, x$y, pch = 21, bg = "white", lwd = 2)
    }
}
