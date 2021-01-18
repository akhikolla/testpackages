.corLevLowRankR <- function(par,
                            nlevels,
                            rank,
                            lowerSQRT = FALSE,
                            compGrad = FALSE,
                            cov = 0) {
    
    
    cov <- as.integer(cov)

    ## for simplicity
    m <- nlevels
    r <- rank
    
    np <- length(par)

    if (cov == 0) {
        npCor <- np
    } else if (cov == 1) {
        npCor <- np - 1
        sigma2 <- par[np]
        par <- par[-np]
        np <- np - 1L
    } else if (cov == 2) {
        npCor <- np - nlevels
        sigma2 <- par[(npCor + 1L):np]
        par <- par[1L:npCor ]
        np <- npCor
    } 
    
    if (npCor != (r - 1L) * (m - r / 2)) {
        stop("bad number of parameters")
    }
   
    Thetas <- matrix(0, nrow = r - 1, ncol = m)
    Thetas[1L, 1L] <- 0

    Thetas[upper.tri(Thetas)] <- par
    Thetas <- t(Thetas)
    
    L <- matrix(0, nrow = m, ncol = r)
    L[1L, 1L] <- 1.0
    
    if (r == 1L) {
        for (i in 2L:m) {
            L[i, 1L] <- cos(Thetas[i, j])
        }
    } else {
        ## case min(i, j) == i
        for (i in 2L:r) {
            prod_sin <- 1
            for (j in 1L:(i - 1L)) {
                L[i, j] <- cos(Thetas[i, j]) * prod_sin
                prod_sin <- prod_sin * sin(Thetas[i, j])
            }
            L[i, i] <- prod_sin
        }
        ## case min(i, j) == r
        if (r < m) {
            for (i in (r + 1L):m) {
                prod_sin <- 1
                for (j in 1L:(r - 1L)) {
                    L[i, j] <- cos(Thetas[i, j]) * prod_sin
                    prod_sin <- prod_sin * sin(Thetas[i, j])
                }
                L[i, r] <- prod_sin
            }
        }
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
##' Compute the correlation matrix for a low-rank structure.
##' 
##' The correlation matrix with size \eqn{m} is the general symmetric
##' correlation matrix with rank \eqn{\leq r}{<= r} where \eqn{r} is
##' given, as described by Rapisarda et al. It depends on \eqn{(r - 1)
##' \times (m - r / 2) / 2}{(r - 1) * (m - r / 2)} parameters
##' \eqn{\theta_{ij}}{theta[i, j]} where the indices \eqn{i} and
##' \eqn{j} are such that \eqn{1 \leq j < i}{1 <= j < i} for \eqn{i
##' \leq r}{i <= r} or such that \eqn{1 \leq j < r}{1 <= j < r} for
##' \eqn{r < i \leq n}{r < i <= m}.  The parameters
##' \eqn{\theta_{ij}}{theta[i, j]} are angles and are to be taken to
##' be in \eqn{[0, 2\pi)}{[0, 2*pi)} if \eqn{j = 1}{j = 1} and in
##' \eqn{[0, \pi)}{[0, pi)} otherwise.
##' 
##' @title Correlation Matrix for a Low-Rank Structure
##'
##' @param par A numeric vector with length \code{npCor + npVar} where
##' \code{npCor = (rank - 1) * (nlevels - rank / 2)} is the number
##' of correlation parameters, and \code{npVar} is the number of
##' variance parameters, which depends on the value of \code{cov}. The
##' value of \code{npVar} is \code{0}, \code{1} or \code{nlevels}
##' corresponding to the values of \code{cov}: \code{0}, \code{1} and
##' \code{2}. The correlation parameters are assumed to be located at
##' the head of \code{par} i.e. at indices \code{1} to
##' \code{npCor}. The variance parameter(s) are assumed to be at the
##' tail, i.e. at indices \code{npCor +1 } to \code{npCor + npVar}.
##'
##' @param nlevels Number of levels \eqn{m}. 
##'
##' @param rank The rank, which must be \eqn{>1} and \eqn{< nlevels}.
##' 
##' @param levels Character representing the levels.
##'
##' @param lowerSQRT Logical. When \code{TRUE} a lower-triangular root
##' \eqn{\mathbf{L}}{L} of the correlation or covariance matrix
##' \eqn{\mathbf{C}}{C} is returned instead of the correlation matrix.
##' Note that this matrix can have negative diagonal elements hence is
##' not a (pivoted) Cholesky root.
##'
##' @param compGrad Logical. Should the gradient be computed? This is
##' only possible for the C implementation.
##'
##' @param cov Integer \code{0}, \code{1} or \code{2}. If \code{cov}
##' is \code{0}, the matrix is a \emph{correlation} matrix (or its
##' root). If \code{cov} is \code{1} or \code{2}, the matrix is a
##' \emph{covariance} (or its root) with constant variance vector for
##' \code{code = 1} and with arbitrary variance for \code{code =
##' 2}. The variance parameters \code{par} are located at the tail of
##' the \code{par} vector, so at locations \code{npCor + 1} to
##' \code{npCor + nlevels} when \code{code = 2} where \code{npCor} is
##' the number of correlation parameters.
##' 
##' @param impl A character telling which of the C and R implementations
##' should be chosen. Not used for now.
##' 
##' @return A correlation matrix (or its root) with the
##' optional \code{gradient} attribute. 
##'
##' @author Yves Deville
##'
##' @note Here the parameters \eqn{\theta_{ij}}{theta[i, j]} are used
##' \emph{in row order} rather than in the column order. This order
##' simplifies the computation of the gradients.
##'
##' @references
##'
##' Francesco Rapisarda, Damanio Brigo, Fabio Mercurio (2007).
##' "Parameterizing Correlations a Geometric Interpretation". \emph{IMA
##' Journal of Management Mathematics}, \bold{18}(1): 55-73.
##'
##' Igor Grubisic, Raoul Pietersz (2007). "Efficient Rank Reduction of
##' Correlation Matrices". \emph{Linear Algebra and its Applications},
##' \bold{422}: 629-653.
##' 
##' @seealso The \code{\link{corLevSymm}} function for the full-rank
##' case.
##' 
##' @examples
corLevLowRank <- function(par,
                          nlevels,
                          rank,
                          levels,
                          lowerSQRT = FALSE,
                          compGrad = TRUE,
                          cov = 0,
                          impl = c("C", "R")) {
    
    if ((rank < 2) || (rank >= nlevels)) {
        stop("'rank' must be >= 2 and < nlevels")     
    }
   
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
        
        res <- .corLevLowRankR(par = par,
                               nlevels = nlevels,
                               rank = rank,
                               lowerSQRT = lowerSQRT,
                               compGrad = compGrad,
                               cov = cov)
        return(res)
    }
   
    np <- length(par)
    npCor <- as.integer((rank - 1) * (nlevels - rank / 2))
    npVar <- c(0L, 1L, nlevels)[cov + 1L]
    
    ## ========================================================================
    ## Manage the choice in 'cov'
    ## ========================================================================

    if (!np) {
        parCor <- rep(pi / 2, npCor)
        names(parCor) <- parNamesSymm(nlevels)
        if (cov == 0) {
            par <- parCor
        } else if (cov == 1) {
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
        corMat <- .Call(corLev_LowRank,
                        parCor,
                        as.integer(nlevels),
                        as.integer(rank),
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
##' low-rank correlation.
##'
##' @title Qualitative Correlation or Covariance Kernel with one Input
##' and Low-Rank Correlation.
##' 
##' @param factor A factor with the wanted levels for the covariance
##' kernel object.
##'
##' @param rank The wanted rank, which must be \eqn{\geq 2}{>= 2} and
##' \eqn{< m} where \eqn{m} is the number of levels.
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
##' @seealso The \code{\link{q1Symm}} function to create a kernel
##' object for the full-rank case.
##' 
##' @examples
q1LowRank <- function(factor,
                      rank = 2L,
                      input = "x",
                      cov = c("corr", "homo", "hete"),
                      intAsChar = TRUE){

    cov <- match.arg(cov)
    cov <- match(cov, c("corr", "homo", "hete")) - 1L
 
    nlev <- nlevels(factor)
    parN <- as.integer((rank - 1L) * (nlev - rank / 2))
    
    inputNames <- input
    lev <- list()
    lev[[inputNames]] <- levels(factor)

    ## find the parameter names: triangular part and rectangular part.            
    kernParNames <- parNamesSymm(rank)

    pnm <- paste("theta",
                 as.vector(t(outer((rank + 1L):nlev, 1L:(rank - 1L),
                                   paste, sep = "_"))),
                 sep = "_")
    
    kernParNames <- c(kernParNames, pnm)

    ## The upper bound on a parameters must be 'pi' except for the
    ## angles in the first column and in rows 2 to 'r"
    
    ii <- 2L:rank
    ii <- ii * (ii - 1) / 2 
    parLower <- rep(0, parN)
    parUpper <- rep(pi, parN)
    parUpper[ii] <- 2 * pi
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
        corLevLowRank(nlevels = nlev,
                      levels = lev[[1]],
                      rank = rank,
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
        label = sprintf("Low rank, rank = %d", rank),
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
