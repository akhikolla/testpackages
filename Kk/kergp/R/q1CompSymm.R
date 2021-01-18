.corLevCompSymmR <- function(par,
                             nlevels,
                             levels,
                             lowerSQRT = FALSE,
                             compGrad = TRUE,
                             cov = 0) {
    
    np <- length(par)
    if (missing(levels)) levels <- 1L:nlevels
    cov <- as.integer(cov)
    
    if (cov) {
        if (np != 2L) {
            stop("parameter vector 'par' with bad length: must be 2")
        }
        sigma2 <- par[2L]
        sigma <- sqrt(sigma2)
    } else {
        if (np != 1L) {
            stop("parameter vector 'par' with bad length: must be 1")
        }
    }
    parNames <- c("rho", "sigma2")[1L:np]
    
    rho <- par[1L]

    if (!lowerSQRT) {
        
        corMat <- array(rho, dim = c(nlevels, nlevels))
        diag(corMat) <- 1.0
        
        if (compGrad) {
            g <- array(1.0, dim = c(nlevels, nlevels))
            diag(g) <- 0.0

            if (cov) {
                g <- c(sigma2 * g, as.vector(corMat))
            }
            dim(g) <- c(nlevels, nlevels, np)
            dimnames(g) <- list(NULL, NULL, parNames)
        } 

        if (cov) corMat <- sigma2 * corMat
        if (compGrad) attr(corMat, "gradient") <- g
        
        
        return(corMat)
    }

    L <- array(0, dim = c(nlevels, nlevels))
    
    a <- rep(NA, nlevels)
    b <- rep(NA, nlevels)
    
    a[1L] <- L[1L, 1L] <- 1.0
    S2 <- 0.0

    if (!compGrad) {
    
        for (j in 1:(nlevels - 1L)) {
            b[j] <- (rho - S2) / a[j]
            S2 <- S2 + b[j]^2
            j1 <- j + 1L
            L[j1, j1] <- a[j1] <- sqrt(1 - S2)
            L[j1:nlevels, j] <- b[j]
        }

        if (cov) L <- sigma * L

        return(L)
        
    } else {
        
        da <- rep(0, nlevels)
        db <- rep(0, nlevels)
        dS2 <- 0.0
        dL <-  array(0, dim = c(nlevels, nlevels, np),
                     dimnames = list(levels, levels, parNames))
        
        for (j in 1:(nlevels - 1L)) {

            b[j] <- (rho - S2) / a[j]
            db[j] <- (1 - dS2 - b[j] * da[j]) / a[j]
            
            S2 <- S2 + b[j]^2
            dS2 <- dS2 + 2 * b[j] * db[j]

            j1 <- j + 1L
            a[j1] <- sqrt(1 - S2)
            L[j1, j1] <- a[j1]
            
            da[j1] <- - dS2 / a[j1] / 2
            dL[j1, j1, 1L] <- da[j1]
            
            L[j1:nlevels, j] <- b[j]
            dL[j1:nlevels, j, 1L] <- db[j]

        }
        
        if (cov) {
            L <- sigma * L
            dL[ , , 1L] <- sigma * dL[ , , 1L]
            dL[ , , 2L] <- L / 2 / sigma
        }
        
        attr(L, "gradient") <- dL
    }
    
    L

}

##' Compute the correlation matrix for a the compound symmetry structure.
##'
##' @title Correlation Matrix for the Compound Symmetry Structure
##'
##' @param par Numeric vector of length \code{1} if \code{cov} is
##' \code{TRUE} or with length \code{2} else. The first element is the
##' correlation coefficient and the second one (when it exists) is the
##' variance.
##' 
##' @param nlevels Number of levels. 
##'
##' @param levels Character representing the levels.
##'
##' @param lowerSQRT Logical. When \code{TRUE} the (lower) Cholesky
##' root \eqn{\mathbf{L}}{L} of the correlation matrix
##' \eqn{\mathbf{C}}{C} is returned instead of the correlation matrix.
##'
##' @param compGrad Logical. Should the gradient be computed? 
##'
##' @param cov Logical. If \code{TRUE} the matrix is a covariance
##' matrix (or its Cholesky root) rather than a correlation matrix and
##' the last element in \code{par} is the variance.
##' 
##' @param impl A character telling which of the C and R implementations
##' should be chosen.
##'
##' @return A correlation matrix (or its Cholesky root) with the
##' optional \code{gradient} attribute.
##'
##' @author Yves Deville
##'
##' @note When \code{lowerSQRT} is \code{FALSE}, the implementation
##' used is always in R because no gain would then result from an
##' implementation in C.
##' 
##' @examples
##' checkGrad <- TRUE
##' lowerSQRT <- FALSE
##' nlevels <- 12
##' set.seed(1234)
##' par <- runif(1L, min = 0, max = pi)
##' 
##' ##============================================================================
##' ## Compare R and C implementations for 'lowerSQRT = TRUE'
##' ##============================================================================
##' tR <- system.time(TR <- corLevCompSymm(nlevels = nlevels, par = par,
##'                                        lowerSQRT = lowerSQRT, impl = "R"))
##' tC <- system.time(T <- corLevCompSymm(nlevels = nlevels, par = par,
##'                                       lowerSQRT = lowerSQRT))
##' tC2 <- system.time(T2 <- corLevCompSymm(nlevels = nlevels, par = par,
##'                                         lowerSQRT = lowerSQRT, compGrad = FALSE))
##' ## time
##' rbind(R = tR, C = tC, C2 = tC2)
##'
##' ## results
##' max(abs(T - TR))
##' max(abs(T2 - TR))
##'
##' ##===========================================================================
##' ## Compare the gradients
##' ##===========================================================================
##'
##' if (checkGrad) {
##'
##'     library(numDeriv)
##'
##'     ##=======================
##'     ## lower SQRT case only
##'     ##========================
##'     JR <- jacobian(fun = corLevCompSymm, x = par, nlevels = nlevels,
##'                    lowerSQRT = lowerSQRT, impl = "R", method = "complex")
##'     J <- attr(T, "gradient")
##'
##'     ## redim and compare.
##'     dim(JR) <- dim(J)
##'     max(abs(J - JR))
##'     nG <- length(JR)
##'     plot(1:nG, as.vector(JR), type = "p", pch = 21, col = "SpringGreen3",
##'          cex = 0.8, ylim = range(J, JR),
##'          main = paste("gradient check, lowerSQRT =", lowerSQRT))
##'     points(x = 1:nG, y = as.vector(J), pch = 16, cex = 0.6, col = "orangered")
##' }
##'
corLevCompSymm <- function(par,
                           nlevels,
                           levels,
                           lowerSQRT = FALSE,
                           compGrad = TRUE,
                           cov = FALSE,
                           impl = c("C", "R")) {
  
    impl <- match.arg(impl)
    if (missing(levels)) levels <- 1L:nlevels
    
    np <- length(par)
    
    if (cov) {
        if (np != 2L) {
            stop("parameter vector 'par' with bad length: must be 2")
        }
        sigma2 <- par[2L]
        sigma <- sqrt(sigma2)
    } else {
        if (np != 1L) {
            stop("parameter vector 'par' with bad length: must be 1")
        }
    }

    rho <- par[1L]
    parNames <- c("rho", "sigma2")[1L:np]
    
    if (!lowerSQRT) {
        
        corMat <- array(rho, dim = c(nlevels, nlevels))
        diag(corMat) <- 1.0
        rownames(corMat) <- colnames(corMat) <- levels
        
        if (compGrad) {
            g <- array(1.0, dim = c(nlevels, nlevels))
            diag(g) <- 0.0

            if (cov) {
                g <- c(sigma2 * g, as.vector(corMat))
            }
            dim(g) <- c(nlevels, nlevels, np)
            dimnames(g) <- list(levels, levels, parNames)
        }

        if (cov) corMat <- sigma2 * corMat
        if (compGrad)  attr(corMat, "gradient") <- g

        rownames(corMat) <- colnames(corMat) <- levels
        
        return(corMat)
        
    } else {

        if (impl == "R")  {
            corMat <- .corLevCompSymmR(par, nlevels, lowerSQRT = TRUE,
                                       compGrad = compGrad, cov = cov)
            rownames(corMat) <- colnames(corMat) <- levels
            return(corMat)
        } else {
            corMat <- .Call(corLev_CompSymm,
                            as.double(rho),
                            as.integer(nlevels), 
                            as.integer(lowerSQRT),
                            as.integer(compGrad))

            rownames(corMat) <- colnames(corMat) <- levels
 
            ## redim. Could certainly be done via .Call, but not so clear in
            ## the doc.
            
            if (compGrad) {

                g <- attr(corMat, "gradient")
                
                if (cov) {
                    g <- c(sigma * g, as.vector(corMat) / 2 / sigma)
                }
                
                dim(g) <- c(nlevels, nlevels, np)
                dimnames(g) <- list(levels, levels, parNames)
                attr(corMat, "gradient") <- g
              
            }
            
            if (cov) corMat <- sigma * corMat
            if (compGrad)  attr(corMat, "gradient") <- g

            corMat
            
        }

    }
    
}

##==============================================================================
##' Qualitative correlation or covariance kernel with one input and
##' compound symmetric correlation.
##'
##' @title Qualitative Correlation or Covariance Kernel with one Input
##' and Compound Symmetric Correlation.
##' 
##' @param factor A factor with the wanted levels for the covariance
##' kernel object.
##' 
##' @param input Name of (qualitative) input for the kernel.
##'
##' @param cov Logical. If \code{TRUE}, the kernel will be a covariance
##' kernel, else it will be a correlation kernel.
##' 
##' @return An object with class \code{"covQual"} with \code{d = 1}
##' qualitative input.
##'
##' @note Correlation kernels are needed in tensor products because
##' the tensor product of two covariance kernels each with unknown
##' variance would not be identifiable.
##'
##' @seealso The \code{\link{corLevCompSymm}} function used to compute
##' the correlation matrix and its gradients w.r.t. the correlation
##' parameters.
##' 
##' @examples
##' School <- factor(1L:3L, labels = c("Bad", "Mean" , "Good"))
##' myCor <- q1CompSymm(School, input = "School")
##' coef(myCor) <- 0.26
##'
##' ## Use a data.frame with a factor
##' set.seed(246)
##' newSchool <- factor(sample(1L:3L, size = 20, replace = TRUE),
##'                     labels = c("Bad", "Mean" , "Good"))
##' C1 <- covMat(myCor, X = data.frame(School = newSchool),
##'              compGrad = FALSE, lowerSQRT = FALSE)
##'
##' ## or use a matrix with integer content
##' C2 <- covMat(myCor, X = cbind(School = as.integer(newSchool)),
##'              compGrad = FALSE, lowerSQRT = FALSE)
##' max(abs(C1 - C2))
q1CompSymm <- function(factor,
                       input = "x",
                       cov = c("corr", "homo"),
                       intAsChar = TRUE) {

    cov <- match.arg(cov)
    cov <- match(cov, c("corr", "homo")) - 1L
  
    nlev <- nlevels(factor)
    parN <- 1L
    inputNames <- input
    lev <- list()
    lev[[inputNames]] <- levels(factor)
 
    kernParNames <- "rho"
    parLower <- rep(-1 / (nlev - 1), parN)
    parUpper <- rep(1.0, parN)
    par <- rep(0.5, parN)

    if (cov) {
        kernParNames <- c(kernParNames, "sigma2")
        parLower <- c(parLower, 0.0)
        parUpper <- c(parUpper, Inf)
        par <- c(par, 1.0)
        parN <- parN + 1L
    }

    names(parLower) <- names(parUpper) <- names(par) <- kernParNames
    
    thisCovAll <- function(par, lowerSQRT = FALSE, compGrad = FALSE) {
        corLevCompSymm(nlevels = nlev,
                       levels = lev[[1]],
                       par = par,
                       lowerSQRT = lowerSQRT,
                       compGrad = compGrad,
                       impl = "C",
                       cov = cov)
    }

    ## the slot 'covLevMat' was added on 2017-10-11 to avoid the
    ## recomputation of the matrix at each 'show' or 'covMat' with
    ## 'X' missing
    
    covLevMat <- thisCovAll(par, lowerSQRT = FALSE, compGrad = FALSE)
    rownames(covLevMat) <- colnames(covLevMat) <- lev[[1L]]
    
    new("covQual",            
        covLevels = thisCovAll,
        covLevMat = covLevMat,
        hasGrad = TRUE,
        acceptLowerSQRT = TRUE,
        label = "Compound Symm.",
        d = 1L,
        inputNames = input,
        nlevels = nlev,
        levels = lev,
        parLower = parLower,
        parUpper = parUpper,
        par = par,
        parN = parN,
        kernParNames = kernParNames,
        ordered = FALSE,
        intAsChar = intAsChar)
    
}
