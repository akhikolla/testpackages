## R implementation for checks (not exported).

.corLevSymmR <- function(par, nlevels,
                         lowerSQRT = TRUE,
                         compGrad = FALSE,
                         cov = 0) {
    
    cov <- as.integer(cov)
    

    np <- length(par)
    if (cov == 1) {
        sigma2 <- par[np]
        par <- par[-np]
        np <- np - 1L
    } else if (cov == 2) {
        npCor <- nlevels * (nlevels - 1) / 2
        sigma2 <- par[(npCor + 1):np]
        par <- par[1L:npCor]
        np <- npCor
    }
    m <- sqrt(2 * np + 1 / 4) + 1 / 2
    
    L <- Thetas <- matrix(0, m, m)
    Thetas[1L, 1L] <- 0; L[1L, 1L] <- 1; 

    Thetas[upper.tri(Thetas)] <- par
    Thetas <- t(Thetas)
    
    for (i in 2L:m) {
        prod_sin <- 1
        for (j in 1L:(i - 1L)) {
            L[i, j] <- cos(Thetas[i, j]) * prod_sin
            prod_sin <- prod_sin * sin(Thetas[i, j])
        }
        L[i, i] <- prod_sin
    }

    if (!lowerSQRT) {
        L  <- tcrossprod(L)
        if (cov) {
            L <- L * tcrossprod(sqrt(sigma2))
        }
    } else if (cov) {
        L <- diag(sqrt(sigma2), nrow = m) %*% L
    }     
    
    L

}

##==============================================================================
##' Compute the correlation matrix for a general symmetric correlation
##' structure.
##' 
##' The correlation matrix with dimension \eqn{n} is the
##' \emph{general symmetric correlation matrix} as described by
##' Pinheiro and Bates and implemented in the \bold{nlme} package. It
##' depends on \eqn{n \times (n - 1) / 2}{n * (n - 1) / 2} parameters
##' \eqn{\theta_{ij}}{theta[i, j]} where the indices \eqn{i} and
##' \eqn{j} are such that \eqn{1 \leq j < i \leq n}{1 <= j < i <= n}.
##' The parameters \eqn{\theta_{ij}}{theta[i, j]} are angles and are to
##' be taken to be in \eqn{[0, \pi)}{[0, pi)}  for a one-to-one
##' parameterisation.
##' 
##' @title Correlation Matrix for a General Symmetric Correlation
##' Structure
##'
##' @param par A numeric vector with length \code{npCor + npVar} where
##' \code{npCor = nlevels * (nlevels - 1) / 2} is the number
##' of correlation parameters, and \code{npVar} is the number of
##' variance parameters, which depends on the value of \code{cov}. The
##' value of \code{npVar} is \code{0}, \code{1} or \code{nlevels}
##' corresponding to the values of \code{cov}: \code{0}, \code{1} and
##' \code{2}. The correlation parameters are assumed to be located at
##' the head of \code{par} i.e. at indices \code{1} to
##' \code{npCor}. The variance parameter(s) are assumed to be at the
##' tail, i.e. at indices \code{npCor +1 } to \code{npCor + npVar}.
##'
##' @param nlevels Number of levels. 
##'
##' @param levels Character representing the levels.
##'
##' @param lowerSQRT Logical. When \code{TRUE} the (lower) Cholesky
##' root \eqn{\mathbf{L}}{L} of the correlation or covariance matrix
##' \eqn{\mathbf{C}}{C} is returned instead of the correlation matrix.
##'
##' @param compGrad Logical. Should the gradient be computed? This is
##' only possible for the C implementation.
##'
##' @param cov Integer \code{0}, \code{1} or \code{2}. If \code{cov}
##' is \code{0}, the matrix is a \emph{correlation} matrix (or its
##' Cholesky root). If \code{cov} is \code{1} or \code{2}, the matrix
##' is a \emph{covariance} (or its Cholesky root) with constant
##' variance vector for \code{code = 1} and with arbitrary variance
##' for \code{code = 2}. The variance parameters \code{par} are
##' located at the tail of the \code{par} vector, so at locations
##' \code{npCor + 1} to \code{npCor + nlevels} when \code{code = 2}
##' where \code{npCor} is the number of correlation parameters,
##' i.e. \code{nlevels * (nlevels - 1) / 2}.
##' 
##' @param impl A character telling which of the C and R implementations
##' should be chosen.
##' 
##' @return A correlation matrix (or its Cholesky root) with the
##' optional \code{gradient} attribute. 
##'
##' @author Yves Deville
##'
##' @note Here the parameters \eqn{\theta_{ij}}{theta[i, j]} are used
##' \emph{in row order} rather than in the column order as in the
##' reference or in the \bold{nlme} package. This order simplifies the
##' computation of the gradients.
##'
##' @references
##' 
##' Jose C. Pinheiro and Douglas M. Bates (1996). "Unconstrained
##' Parameterizations for Variance-Covariance matrices". \emph{Statistics and
##' Computing}, 6(3) pp. 289-296.
##'
##' Jose C. Pinheiro and Douglas M. Bates (2000) \emph{Mixed-Effects
##' Models in S and S-PLUS}, Springer.
##' 
##' @seealso The \code{\link{corSymm}} correlation structure in the \bold{nlme}
##' package.
##' 
##' @examples
##' checkGrad <- TRUE
##' nlevels <- 12
##' npar <- nlevels * (nlevels - 1) / 2
##' par <- runif(npar, min = 0, max = pi)
##' ##============================================================================
##' ## Compare R and C implementations for 'lowerSQRT = TRUE'
##' ##============================================================================
##' tR <- system.time(TR <- corLevSymm(nlevels = nlevels,
##'                                    par = par, lowerSQRT = TRUE, impl = "R"))
##' tC <- system.time(T <- corLevSymm(nlevels = nlevels, par = par,
##'                                   lowerSQRT = TRUE))
##' tC2 <- system.time(T2 <- corLevSymm(nlevels = nlevels, par = par, lowerSQRT = TRUE,
##'                                     compGrad = FALSE))
##' ## time
##' rbind(R = tR, C = tC, C2 = tC2)
##'
##' ## results
##' max(abs(T - TR))
##' max(abs(T2 - TR))
##'
##' ##============================================================================
##' ## Compare R and C implementations for 'lowerSQRT = FALSE'
##' ##============================================================================
##' tR <- system.time(TRF <- corLevSymm(nlevels = nlevels, par = par,
##'                                     lowerSQRT = FALSE, impl = "R"))
##' tC <- system.time(TCF <- corLevSymm(nlevels = nlevels, par = par,
##'                                     compGrad = FALSE, lowerSQRT = FALSE))
##' tC2 <- system.time(TCF2 <- corLevSymm(nlevels = nlevels, par = par,
##'                                       compGrad = TRUE, lowerSQRT = FALSE))
##' rbind(R = tR, C = tC, C2 = tC2)
##' max(abs(TCF - TRF))
##' max(abs(TCF2 - TRF))
##'
##' ##===========================================================================
##' ## Compare the gradients
##' ##============================================================================
##'
##' if (checkGrad) {
##'
##'     library(numDeriv)
##'
##'     ##==================
##'     ## lower SQRT case
##'     ##==================
##'     JR <- jacobian(fun = corLevSymm, x = par, nlevels = nlevels, lowerSQRT = TRUE,
##'                    method = "complex", impl = "R")
##'     J <- attr(T, "gradient")
##'
##'     ## redim and compare.
##'     dim(JR) <- dim(J)
##'     max(abs(J - JR))
##'     nG <- length(JR)
##'     plot(1:nG, as.vector(JR), type = "p", pch = 21, col = "SpringGreen3",
##'          cex = 0.8, ylim = range(J, JR),
##'          main = "gradient check, lowerSQRT = TRUE")
##'     points(x = 1:nG, y = as.vector(J), pch = 16, cex = 0.6, col = "orangered")
##'
##'     ##==================
##'     ## Symmetric case
##'     ##==================
##'     JR <- jacobian(fun = corLevSymm, x = par, nlevels = nlevels,
##'                    lowerSQRT = FALSE, impl = "R", method = "complex")
##'     J <- attr(TCF2, "gradient")
##'
##'     ## redim and compare.
##'     dim(JR) <- dim(J)
##'     max(abs(J - JR))
##'     nG <- length(JR)
##'     plot(1:nG, as.vector(JR), type = "p", pch = 21, col = "SpringGreen3",
##'          cex = 0.8,
##'          ylim = range(J, JR),
##'          main = "gradient check, lowerSQRT = FALSE")
##'     points(x = 1:nG, y = as.vector(J), pch = 16, cex = 0.6, col = "orangered")
##' }
corLevSymm <- function(par,
                       nlevels,
                       levels,
                       lowerSQRT = FALSE,
                       compGrad = TRUE,
                       cov = 0,
                       impl = c("C", "R")) {

    cov <- as.integer(cov)
    if (!cov %in% c(0L, 1L, 2L)) {
        stop("'cov' must be 0, 1 or 2")
    }
    

    if (missing(levels)) levels <- 1L:nlevels
    else nlevels <- length(levels)
    
    impl <- match.arg(impl)
    
    if (impl == "R") {

        if (length(par) == 0) {
            stop("the R inplementation can be used only ",
                 "when length(par) > 0")
        }
        
        res <- .corLevSymmR(par, nlevels, lowerSQRT = lowerSQRT, compGrad = compGrad,
                            cov = cov)
        return(res)
    }
   
    np <- length(par)
    npCor <- nlevels * (nlevels - 1L) / 2
    npVar <- c(0, 1, nlevels)[cov + 1L]
    
    ## ========================================================================
    ## Manage the choice in 'cov'
    ## ========================================================================

    if (!np) {
        parCor <- rep(pi / 2, npCor)
        names(parCor) <- parNamesSymm(nlevels)
        if (cov == 0) {
            par <- parCor
        } else  if (cov == 1) {
            sigma2 <- 1.0
            sigma <- 1.0
            par <- c(par, "sigma2" = sigma2)
        } else if (cov == 2) {
            parVar <- rep(1.0, nlevels)
            names(parVar) <- paste("sigma2", 1L:nlevels, sep = "_")
            par <- c(par, parVar)
        }
        np <- length(par)
    } else {
        
        if (np != npCor + npVar) {
            stop("'par' must be of length ", npCor, " + ", npVar,
                 " but is ", np)
        }
        if (cov) {
            sigma2 <- par[(npCor + 1):np]
            sigma <- sqrt(sigma2)
            parCor <- par[1L:npCor]
        } else {
            parCor <- par
        }
    }

    if (npCor) {
        corMat <- .Call(corLev_Symm,
                        parCor,
                        as.integer(nlevels), 
                        as.integer(lowerSQRT),
                        as.integer(compGrad))
    } else {
        corMat <- matrix(1.0, nrow = 1L, ncol = 1L)
        if (cov == 0) {
            rownames(corMat) <- colnames(corMat) <- levels
            return(corMat)
        }
    }
        
    ## ========================================================================
    ## from now on, 'covMat' will be the COVARIANCE matrix or its
    ## lowqer sqrt. It is simpler to keep both matrices alive to cope
    ## with the gradient.
    ## ========================================================================
    
    if (cov == 0) {
        covMat <- corMat
    } else if (cov == 1) {
        if (lowerSQRT) covMat <- sigma * corMat
        else covMat <- sigma2 * corMat
    } else if (cov == 2) {
        if (lowerSQRT) {
            covMat <- sweep(corMat, MARGIN = 1, STATS = sigma, FUN = "*")
        } else {
            sigsig <- tcrossprod(sigma)
            covMat <- sigsig * corMat
        }
    }

    rownames(covMat) <- colnames(covMat) <- levels
    
    ## ========================================================================
    ## redim. Could certainly be done via .Call, but not so clear in
    ## the doc.
    ## ========================================================================
    
    if (compGrad) {
        g <- attr(corMat, "gradient")
                            
        if (cov == 1) {
            if (lowerSQRT) {
                ## XXXX there was an error here: a sqrt was missing!
                g <- c(sigma * g, as.vector(corMat) / 2 / sigma)
            } else {
                g <- c(sigma2 * g, as.vector(corMat))
            }
            dim(g) <- c(nlevels, nlevels, np)
            dimnames(g) <- list(levels, levels, names(par))
        } else if (cov == 2) {

            gBak <- g
            dim(gBak) <- c(nlevels, nlevels, npCor)
            g <- array(0.0, dim = c(nlevels, nlevels, np))
            dimnames(g) <- list(levels, levels, names(par))
            
            if (lowerSQRT) {
                ## this is equivalent to left mulitplying each
                ## gBak[ , , i] by diag(sigma, nrow = nlevels)
                g[ , , 1L:npCor] <- sweep(gBak, MARGIN = 1L, STATS = sigma,
                                          FUN = "*") 
                for (i in 1L:npVar) {
                    g[i, , i + npCor] <- covMat[i, , drop = FALSE] / 2 / sigma2[i]
                }
            } else {
                g[ , , 1L:npCor] <- sweep(gBak, MARGIN = c(1L, 2L),
                                          STATS = sigsig, FUN = "*")
                for (i in 1L:npVar) {
                    prov <- covMat[i, , drop = FALSE] / 2 / sigma2[i]
                    g[ , i, i + npCor] <- prov
                    g[i, , i + npCor] <-  g[i, , i + npCor] + prov
                }
            }
        } else {
            dim(g) <- c(nlevels, nlevels, np)
            dimnames(g) <- list(levels, levels, names(par))
        }
    }
    
    if (compGrad) {
        attr(covMat, "gradient") <- g
    }
    
    covMat
    
}

##==============================================================================
##' Qualitative correlation or covariance kernel with one input and
##' general symmetric correlation.
##'
##' @title Qualitative Correlation or Covariance Kernel with one Input
##' and General Symmetric Correlation.
##' 
##' @param factor A factor with the wanted levels for the covariance
##' kernel object.
##' 
##' @param input Name of (qualitative) input for the kernel.
##'
##' @param cov Integer with value \code{0}, \code{1} or \code{2}. The
##' value \code{0} corresponds to a correlation kernel. The values
##' \code{1} and \code{2} correspond to covariance kernel, with
##' constant variance for \code{cov = 1}, and arbitrary variance
##' vector for \code{cov = 2}.
##' 
##' @return An object with class \code{"covQual"} with \code{d = 1}
##' qualitative input.
##'
##' @note Correlation kernels are needed in tensor products because
##' the tensor product of two covariance kernels each with unknown
##' variance would not be identifiable.
##'
##' @seealso The \code{\link{corLevSymm}} function used to compute the
##' correlation matrix and its gradients w.r.t. the correlation
##' parameters.
##' 
##' @examples
##' School <- factor(1L:3L, labels = c("Bad", "Mean" , "Good"))
##' myCor <- q1Symm(School, input = "School")
##' coef(myCor) <- c(theta_2_1 = pi / 3, theta_3_1 = pi / 4, theta_3_2 = pi / 8)
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
q1Symm <- function(factor, input = "x", cov = c("corr", "homo", "hete"),
                   intAsChar = TRUE){

    cov <- match.arg(cov)
    cov <- match(cov, c("corr", "homo", "hete")) - 1L
 
    nlev <- nlevels(factor)
    
    parN <- as.integer(nlev  * (nlev - 1L) / 2)
    
    inputNames <- input
    lev <- list()
    lev[[inputNames]] <- levels(factor)
 
    kernParNames <- parNamesSymm(nlev)
    parLower <- rep(0, parN)
    parUpper <- rep(pi, parN)
    par <- rep(pi / 2, parN)

    if (cov == 1L) {
        kernParNames <- c(kernParNames, "sigma2")
        parLower <- c(parLower, 0.0)
        parUpper <- c(parUpper, Inf)
        par <- c(par, 1.0)
        parN <- parN + 1L
    } else if (cov == 2L) {
        kernParNames <- c(kernParNames, paste("sigma2", 1L:nlev, sep = "_"))
        parLower <- c(parLower, rep(0.0, nlev))
        parUpper <- c(parUpper, rep(Inf, nlev))
        par <- c(par, rep(1.0, nlev))
        parN <- parN + nlev
    }
    
    if (parN) {
        names(parLower) <- names(parUpper) <- names(par) <- kernParNames
    }

    thisCovAll <- function(par, lowerSQRT = FALSE, compGrad = FALSE) {
        corLevSymm(nlevels = nlev,
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
        label = "General Symm.",
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
