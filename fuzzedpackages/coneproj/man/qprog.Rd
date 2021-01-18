\name{qprog}
\alias{qprog}
\title{Quadratic Programming}
\description{
Given a positive definite \eqn{n} by \eqn{n} matrix \eqn{Q} and a constant vector \eqn{c} in \eqn{R^n}, the object is to find \eqn{\theta} in \eqn{R^n} to minimize \eqn{\theta'Q\theta - 2c'\theta} subject to \eqn{A\theta \ge b}, for an irreducible constraint matrix \eqn{A}. This routine transforms into a cone projection problem for the constrained solution. }
\usage{qprog(q, c, amat, b, face = NULL, msg = TRUE)}
\arguments{
  \item{q}{A \eqn{n} by \eqn{n} positive definite matrix.}
  \item{c}{A vector of length \eqn{n}.}
  \item{amat}{A \eqn{m} by \eqn{n} constraint matrix. The rows of amat must be irreducible.}
  \item{b}{A vector of length \eqn{m}. Its default value is 0.}
    \item{face}{A vector of the positions of edges, which define the initial face for the cone projection. For example, when there are \eqn{m} cone edges, then face is a subset of \eqn{1,\ldots,m}. The default is face = NULL.}
  \item{msg}{A logical flag. If msg is TRUE, then a warning message will be printed when there is a non-convergence problem; otherwise no warning message will be printed. The default is msg = TRUE}
}
\details{
To get the constrained solution to \eqn{\theta'Q\theta - 2c'\theta} subject to \eqn{A\theta \ge b}, this routine makes the Cholesky decomposition of \eqn{Q}. Let \eqn{U'U = Q}, and define \eqn{\phi = U\theta} and \eqn{z = U^{-1}c}, where \eqn{U^{-1}} is the inverse of \eqn{U}.
Then we minimize \eqn{||z - \phi||^2}, subject to \eqn{B\phi \ge 0}, where \eqn{B = AU^{-1}}. It is now a cone projection problem with the constraint cone \eqn{C} of the form \eqn{\{\phi: B\phi \ge 0 \}}. This routine gives the estimation of \eqn{\theta}, which is \eqn{U^{-1}} times the estimation of \eqn{\phi}. 

The routine qprog dynamically loads a C++ subroutine "qprogCpp".
}

\value{
\item{df}{The dimension of the face of the constraint cone on which the projection lands.}
\item{thetahat}{A vector minimizing \eqn{\theta'Q\theta - 2c'\theta}.}
\item{steps}{The number of iterations before the algorithm converges.}
\item{xmat}{The rows of the matrix are the edges of the face of the polar cone on which the residual of the projection onto the constraint cone lands.}
\item{face}{A vector of the positions of edges, which define the face on which the final projection lands on. For example, when there are \eqn{m} cone edges, then face is a subset of \eqn{1,\ldots,m}.}
}
\references{
Goldfarb, D. and A. Idnani (1983) A numerically stable dual method for solving strictly convex quadratic programs.
\emph{Mathematical Programming \bold{27}}, 1--33.

Fraser, D. A. S. and H. Massam (1989) A mixed primal-dual bases algorithm for regression
under inequality constraints application to concave regression. \emph{Scandinavian Journal of Statistics \bold{16}}, 65--74.

Fang,S.-C. and S. Puthenpura (1993) \emph{Linear Optimization and Extensions}. Englewood Cliffs, New Jersey: Prentice Hall.

Silvapulle, M. J. and P. Sen (2005) \emph{Constrained Statistical Inference}. John Wiley and Sons.

Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications
in Statistics \bold{42(5)}}, 1126--1139.

Liao, X. and M. C. Meyer (2014) coneproj: An R package for the primal or dual cone projections with routines for constrained regression. \emph{Journal of Statistical Software \bold{61(12)}}, 1--22.
}
\author{
Mary C. Meyer and Xiyue Liao 
}

\seealso{
 \code{\link{coneA}}
}
\examples{
# load the cubic data set
    data(cubic)

# extract x
    x <- cubic$x

# extract y
    y <- cubic$y

# make the design matrix
    xmat <- cbind(1, x, x^2, x^3)

# make the q matrix
    q <- crossprod(xmat)

# make the c vector
    c <- crossprod(xmat, y)

# make the constraint matrix to constrain the regression to be increasing, nonnegative and convex
    amat <- matrix(0, 4, 4)
    amat[1, 1] <- 1; amat[2, 2] <- 1
    amat[3, 3] <- 1; amat[4, 3] <- 1
    amat[4, 4] <- 6
    b <- rep(0, 4)

# call qprog 
    ans <- qprog(q, c, amat, b)

# get the constrained fit of y
    betahat <- fitted(ans)
    fitc <- crossprod(t(xmat), betahat)

# get the unconstrained fit of y
    fitu <- lm(y ~ x + I(x^2) + I(x^3))

# make a plot to compare fitc and fitu
    par(mar = c(4, 4, 1, 1))
    plot(x, y, cex = .7, xlab = "x", ylab = "y")
    lines(x, fitted(fitu))
    lines(x, fitc, col = 2, lty = 4)
    legend("topleft", bty = "n", c("constr.fit", "unconstr.fit"), lty = c(4, 1), col = c(2, 1))
    title("Qprog Example Plot")
}

\keyword{cone projection routines}
