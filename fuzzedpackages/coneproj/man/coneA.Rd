\name{coneA}
\alias{coneA}
\title{Cone Projection -- Polar Cone}
\description{
This routine implements the hinge algorithm for cone projection to minimize 
\eqn{||y - \theta||^2} over the cone \eqn{C} of the form \eqn{\{\theta: A\theta \ge 0\}}. }
\usage{coneA(y, amat, w = NULL, face = NULL, msg = TRUE)  }
\arguments{
  \item{y}{A vector of length \eqn{n}.}
  \item{amat}{A constraint matrix. The rows of amat must be irreducible. The column number of amat must equal the length of \eqn{y}.}
  \item{w}{An optional nonnegative vector of weights of length \eqn{n}. If w is not given, all weights are taken to equal 1. Otherwise, the minimization of \eqn{(y - \theta)'w(y - \theta)} over \eqn{C} is returned. The default is w = NULL.}
   \item{face}{A vector of the positions of edges, which define the initial face for the cone projection. For example, when there are \eqn{m} cone edges, then face is a subset of \eqn{1,\ldots,m}. The default is face = NULL.}
  \item{msg}{A logical flag. If msg is TRUE, then a warning message will be printed when there is a non-convergence problem; otherwise no warning message will be printed. The default is msg = TRUE}
}

\details{
The routine coneA dynamically loads a C++ subroutine "coneACpp". The rows
of \eqn{- A} are the edges of the polar cone \eqn{\Omega^o}. This routine first projects \eqn{y} onto \eqn{\Omega^o} to get the residual of the projection onto the constraint cone \eqn{C}, and then uses the fact that \eqn{y} is equal to the sum of the projection of \eqn{y} onto \eqn{C} and the projection of \eqn{y} onto \eqn{\Omega^o} to get the estimation of \eqn{\theta}. See references cited in this section for more details about the relationship between polar cone and constraint cone. }

\value{
\item{df}{The dimension of the face of the constraint cone on which the projection lands.}
\item{thetahat}{The projection of  \eqn{y} on the constraint cone.}
\item{steps}{The number of iterations before the algorithm converges.}
\item{xmat}{The rows of the matrix are the edges of the face of the polar cone on which the residual of the projection onto the constraint cone lands.}
\item{face}{A vector of the positions of edges, which define the face on which the final projection lands on. For example, when there are \eqn{m} cone edges, then face is a subset of \eqn{1,\ldots,m}.}
}

\references{
Meyer, M. C. (1999) An extension of the mixed primal-dual bases algorithm to the case of more constraints than 
dimensions. \emph{Journal of Statistical Planning and Inference \bold{81}}, 13--31.

Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications
in Statistics \bold{42(5)}}, 1126--1139.  

Liao, X. and M. C. Meyer (2014) coneproj: An R package for the primal or dual cone projections with routines for constrained regression. \emph{Journal of Statistical Software \bold{61(12)}}, 1--22.
}

\author{
Mary C. Meyer and Xiyue Liao  }

\seealso{
\code{\link{coneB}}, \code{\link{constreg}}, \code{\link{qprog}}  }

\examples{
# generate y
    set.seed(123)
    n <- 50
    x <- seq(-2, 2, length = 50)
    y <- - x^2 + rnorm(n)

# create the constraint matrix to make the first half of y monotonically increasing
# and the second half of y monotonically decreasing
    amat <- matrix(0, n - 1, n)
    for(i in 1:(n/2 - 1)){
       amat[i, i] <- -1; amat[i, i + 1] <- 1
    }
    for(i in (n/2):(n - 1)){
       amat[i, i] <- 1; amat[i, i + 1] <- -1
    }

# call coneA
    ans1 <- coneA(y, amat)
    ans2 <- coneA(y, amat, w = (1:n)/n)

# make a plot to compare the unweighted fit and the weighted fit
    par(mar = c(4, 4, 1, 1))
    plot(y, cex = .7, ylab = "y")
    lines(fitted(ans1), col = 2, lty = 2)
    lines(fitted(ans2), col = 4, lty = 2)
    legend("topleft", bty = "n", c("unweighted fit", "weighted fit"), col = c(2, 4), lty = c(2, 2))
    title("ConeA Example Plot")  
}

