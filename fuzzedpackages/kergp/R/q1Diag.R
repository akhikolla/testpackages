##==============================================================================
##' Compute the correlation or covariance matrix for a diagonal
##' structure.
##' 
##' @title Correlation or Covariance Matrix for a Diagonal Structure
##'
##' @param par A numeric vector with length \code{npVar} where
##' \code{npVar} is the number of variance parameters, namely
##' \code{0}, \code{1} or \code{nlevels} corresponding to the values
##' of \code{cov}: \code{0}, \code{1} and \code{2}. 
##'
##' @param nlevels Number of levels. 
##'
##' @param levels Character representing the levels.
##'
##' @param lowerSQRT Logical. When \code{TRUE} the (lower) Cholesky
##' root \eqn{\mathbf{L}}{L} of the correlation or covariance matrix
##' \eqn{\mathbf{C}}{C} is returned instead of the correlation matrix.
##'
##' @param compGrad Logical. Should the gradient be computed? 
##'
##' @param cov Integer \code{0}, \code{1} or \code{2}. If \code{cov}
##' is \code{0}, the matrix is a \emph{correlation} matrix (or its
##' Cholesky root) i.e. an identity matrix. If \code{cov} is \code{1}
##' or \code{2}, the matrix is a \emph{covariance} (or its square
##' root) with constant variance vector for \code{code = 1} and with
##' arbitrary variance vector for \code{code = 2}. 
##' 
##' @return A correlation matrix (or its Cholesky root) with the
##' optional \code{gradient} attribute. 
##' 
##' @examples
##' set.seed(123)
##' checkGrad <- TRUE
##' nlevels <- 12
##' sigma2 <- rexp(n = nlevels)
##' T0 <- corLevDiag(nlevels = nlevels, par = sigma2, cov = 2)
##' L0 <- corLevDiag(nlevels = nlevels, par = sigma2, cov = 2,
##'                  lowerSQRT = TRUE)
##' 
corLevDiag <- function(par,
                       nlevels,
                       levels,
                       lowerSQRT = FALSE,
                       compGrad = TRUE,
                       cov = 0) {
    
    if (missing(levels)) {
        if (missing(nlevels)) {
            stop("when 'levels' is missing, 'nlevels' should be given")
        }
        nlevels <- as.integer(nlevels)
        levels <- 1L:nlevels
    } else nlevels <- length(levels)
    
    cov <- as.integer(cov)
    if (!cov %in% c(0L, 1L, 2L)) {
        stop("'cov' must be 0, 1 or 2")
    }

    np <- length(par)
    npVar <- c(0L, 1L, nlevels)[cov + 1L]
    
    ## ========================================================================
    ## Manage the choice in 'cov'. If 'np' is zero we set the parameters to
    ## default values. 
    ## ========================================================================

    if (!np) {
        if (cov == 0) {
            par <- numeric(0)
        } else  if (cov == 1L) {
            par <- c("sigma2" = 1.0)
        } else if (cov == 2) {
            par <- rep(1.0, nlevels)
            names(par) <- paste("sigma2", 1L:nlevels, sep = "_")
        }
        np <- length(par)
    } else if (np != npVar) {
            stop("'par' must be of length ", npVar,
                 " but is ", np)
    }
    
    sigma2 <- par
    sigma <- sqrt(sigma2)

    if (cov == 0L) {
        covMat <- diag(nrow = nlevels)
    } else {
        if (!lowerSQRT) {
            covMat <- diag(x = sigma2, nrow = nlevels)
        } else {
            covMat <- diag(x = sigma, nrow = nlevels)
        }
    }
    
    if (!compGrad) return(covMat)
    
    g <- array(0.0, dim = c(nlevels, nlevels, np),
               dimnames = list(levels, levels, names(par)))

    ## this works for both cases cov = 1 (npar = 1) and cov = 2
    ## (npar = nlevels)
    if (cov) {
        if (!lowerSQRT) {
            for (i in seq_along(par)) {
                g[ , , i] <- diag(nrow = nlevels)
            }
        } else {
            for (i in seq_along(par)) {
                 g[ , , i] <- diag(nrow = nlevels) * 2.0 / sigma[i]
             }
        }
    }

    attr(covMat, "gradient") <- g
    covMat
    
}

##==============================================================================
##' Qualitative correlation or covariance kernel with one input and
##' diagonal structure.
##'
##' @title Qualitative Correlation or Covariance Kernel with one Input
##' and Diagonal Structure
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
##' @note The correlation version obtained with \code{cov = 0} has no
##' parameters.
##'
##' @seealso \code{\link{q1Symm}}, \code{\link{q1CompSymm}} are other
##' covariance structures for one  qualitative input.
##' 
##' @examples
##' School <- factor(1L:3L, labels = c("Bad", "Mean" , "Good"))
##' 
##' ## correlation: no parameter!
##' myCor <- q1Diag(School, input = "School")
##'
##' ## covariance 
##' myCov <- q1Diag(School, input = "School", cov = 2)
##' coef(myCov) <- c(1.1, 2.2, 3.3)
##' 
q1Diag <- function(factor, input = "x", cov = c("corr", "homo", "hete"),
                   intAsChar = TRUE){

    cov <- match.arg(cov)
    cov <- match(cov, c("corr", "homo", "hete")) - 1L
    
    nlev <- nlevels(factor)
    
    inputNames <- input
    lev <- list()
    lev[[inputNames]] <- levels(factor)

    if (cov == 0L) {
        kernParNames <- character(0)
        parLower <- parUpper <- par <- numeric(0)
        parN <- 0L
    } else  if (cov == 1L) {
        kernParNames <- "sigma2"
        parLower <- 0.0
        parUpper <- Inf
        par <- 1.0
        parN <- 1L
    } else if (cov == 2L) {
        kernParNames <- paste("sigma2", 1L:nlev, sep = "_")
        parLower <- rep(0.0, nlev)
        parUpper <- rep(Inf, nlev)
        par <- rep(1.0, nlev)
        parN <- nlev
    }
    
    if (parN) {
        names(parLower) <- names(parUpper) <- names(par) <- kernParNames
    }

    thisCovAll <- function(par, lowerSQRT = FALSE, compGrad = FALSE) {
        corLevDiag(nlevels = nlev,
                   levels = lev[[1L]],
                   par = par,
                   lowerSQRT = lowerSQRT,
                   compGrad = compGrad,
                   cov = cov)
    }

    ## ========================================================================
    ## the slot 'covLevMat' was added on 2017-10-11 to avoid the
    ## recomputation of the matrix at each 'show' or 'covMat' with 'X'
    ## missing
    ## ========================================================================
    
    covLevMat <- thisCovAll(par, lowerSQRT = FALSE, compGrad = FALSE)
    rownames(covLevMat) <- colnames(covLevMat) <- lev[[1L]]
    
    new("covQual",            
        covLevels = thisCovAll,
        covLevMat = covLevMat,
        hasGrad = TRUE,
        acceptLowerSQRT = TRUE,
        label = "Diagonal cov.",
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
