\name{corLevSymm}
\alias{corLevSymm}
\title{

  Correlation Matrix for a General Symmetric Correlation Structure

}
\usage{

corLevSymm(par, nlevels, levels, lowerSQRT = FALSE, compGrad = TRUE,
           cov = 0, impl = c("C", "R"))

}

\arguments{

  \item{par}{

    A numeric vector with length \code{npCor + npVar} where
    \code{npCor = nlevels *} \code{(nlevels - 1) / 2} is the number of
    correlation parameters, and \code{npVar} is the number of variance
    parameters, which depends on the value of \code{cov}. The value of \code{npVar}
    is \code{0}, \code{1} or \code{nlevels} corresponding to the values
    of \code{cov}: \code{0}, \code{1} and \code{2}. The correlation
    parameters are assumed to be located at the head of \code{par}
    i.e. at indices \code{1} to \code{npCor}. The variance parameter(s)
    are assumed to be at the tail, i.e. at indices \code{npCor + 1 } to
    \code{npCor + npVar}.

}
\item{nlevels}{

  Number of levels.

}
\item{levels}{

  Character representing the levels.

}
\item{lowerSQRT}{
  
  Logical. When \code{TRUE} the (lower) Cholesky root
  \eqn{\mathbf{L}}{L} of the correlation or covariance matrix
  \eqn{\mathbf{C}}{C} is returned instead of the correlation matrix.

}

\item{compGrad}{

  Logical. Should the gradient be computed? This is only possible for
  the C implementation.

}
\item{cov}{

  Integer \code{0}, \code{1} or \code{2}. If \code{cov} is \code{0}, the
  matrix is a \emph{correlation} matrix (or its Cholesky root). If
  \code{cov} is \code{1} or \code{2}, the matrix is a \emph{covariance}
  (or its Cholesky root) with constant variance vector for
  \code{code = 1} and with arbitrary variance for \code{code = 2}. The variance
  parameters \code{par} are located at the tail of the \code{par}
  vector, so at locations \code{npCor + 1} to \code{npCor + nlevels}
  when \code{code = 2} where \code{npCor} is the number of correlation
  parameters, i.e. \code{nlevels * (nlevels - 1) / 2}.

}
\item{impl}{

  A character telling which of the C and R implementations should be
  chosen.

}

}
\value{

  A correlation matrix (or its Cholesky root) with the optional
  \code{gradient} attribute.

}
\description{
  
  Compute the correlation matrix for a general symmetric correlation
  structure.
  
}
\details{

  The correlation matrix with dimension \eqn{n} is the \emph{general
    symmetric correlation matrix} as described by Pinheiro and Bates and
  implemented in the \pkg{nlme} package. It depends on \eqn{n \times (n -
    1) / 2}{n * (n - 1) / 2} parameters \eqn{\theta_{ij}}{theta[i, j]} where
  the indices \eqn{i} and \eqn{j} are such that \eqn{1 \leq j < i \leq
    n}{1 <= j < i <= n}.  The parameters \eqn{\theta_{ij}}{theta[i, j]} are
  angles and are to be taken to be in \eqn{[0, \pi)}{[0, pi)} for a
  one-to-one parameterisation.
  
}
\note{

  This function is essentially for internal use and the corresponding
  correlation or covariance kernels are created as \code{covQual}
  objects by using the \code{\link{q1Symm}} creator. 

  The parameters \eqn{\theta_{ij}}{theta[i, j]} are used \emph{in
    row order} rather than in the column order as in the reference or in the
  \pkg{nlme} package. This order simplifies the computation of the
  gradients.

}
\examples{
checkGrad <- TRUE
nlevels <- 12
npar <- nlevels * (nlevels - 1) / 2
par <- runif(npar, min = 0, max = pi)
##============================================================================
## Compare R and C implementations for 'lowerSQRT = TRUE'
##============================================================================
tR <- system.time(TR <- corLevSymm(nlevels = nlevels,
                                   par = par, lowerSQRT = TRUE, impl = "R"))
tC <- system.time(T <- corLevSymm(nlevels = nlevels, par = par,
                                  lowerSQRT = TRUE))
tC2 <- system.time(T2 <- corLevSymm(nlevels = nlevels, par = par,
                                    lowerSQRT = TRUE, compGrad = FALSE))
## time
rbind(R = tR, C = tC, C2 = tC2)

## results
max(abs(T - TR))
max(abs(T2 - TR))

##============================================================================
## Compare R and C implementations for 'lowerSQRT = FALSE'
##============================================================================
tR <- system.time(TRF <- corLevSymm(nlevels = nlevels, par = par,
                                    lowerSQRT = FALSE, impl = "R"))
tC <- system.time(TCF <- corLevSymm(nlevels = nlevels, par = par,
                                    compGrad = FALSE, lowerSQRT = FALSE))
tC2 <- system.time(TCF2 <- corLevSymm(nlevels = nlevels, par = par,
                                      compGrad = TRUE, lowerSQRT = FALSE))
rbind(R = tR, C = tC, C2 = tC2)
max(abs(TCF - TRF))
max(abs(TCF2 - TRF))

##===========================================================================
## Compare the gradients
##===========================================================================

if (checkGrad) {

    library(numDeriv)

    ##==================
    ## lower SQRT case
    ##==================
    JR <- jacobian(fun = corLevSymm, x = par, nlevels = nlevels,
                   lowerSQRT = TRUE, method = "complex", impl = "R")
    J <- attr(T, "gradient")

    ## redim and compare.
    dim(JR) <- dim(J)
    max(abs(J - JR))
    nG <- length(JR)
    plot(1:nG, as.vector(JR), type = "p", pch = 21, col = "SpringGreen3",
         cex = 0.8, ylim = range(J, JR),
         main = "gradient check, lowerSQRT = TRUE")
    points(x = 1:nG, y = as.vector(J), pch = 16, cex = 0.6, col = "orangered")

    ##==================
    ## Symmetric case
    ##==================
    JR <- jacobian(fun = corLevSymm, x = par, nlevels = nlevels,
                   lowerSQRT = FALSE, impl = "R", method = "complex")
    J <- attr(TCF2, "gradient")

    ## redim and compare.
    dim(JR) <- dim(J)
    max(abs(J - JR))
    nG <- length(JR)
    plot(1:nG, as.vector(JR), type = "p", pch = 21, col = "SpringGreen3",
         cex = 0.8,
         ylim = range(J, JR),
         main = "gradient check, lowerSQRT = FALSE")
    points(x = 1:nG, y = as.vector(J), pch = 16, cex = 0.6, col = "orangered")
}
}
\references{
  
Jose C. Pinheiro and Douglas M. Bates
(1996). "Unconstrained Parameterizations for Variance-Covariance matrices".
\emph{Statistics and Computing}, 6(3) pp. 289-296.

Jose C. Pinheiro and Douglas M. Bates (2000) \emph{Mixed-Effects
Models in S and S-PLUS}, Springer.

}
\seealso{
  
  The \code{corSymm} correlation structure in the \pkg{nlme}
  package.

}
%% \author{Yves Deville}
