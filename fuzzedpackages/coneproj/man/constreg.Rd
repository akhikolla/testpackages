\name{constreg}
\alias{constreg}
\title{Constrained Parametric Regression}

\description{
The least-squares regression model \eqn{y = X\beta + \varepsilon} is considered, where the object is to find \eqn{\beta} to minimize
\eqn{||y - X\beta||^2}, subject to \eqn{A\beta \ge 0}.  }

\usage{
constreg(y, xmat, amat, w = NULL, test = FALSE, nloop = 1e+4)
}

\arguments{
  \item{y}{A vector of length \eqn{n}.}
  \item{xmat}{A full column-rank design matrix. The column number of xmat must equal the length of \eqn{\beta}.}
  \item{amat}{A constraint matrix. The rows of amat must be irreducible. The column number of amat must equal the length of \eqn{\beta}.}
 \item{w}{An optional nonnegative vector of weights of length \eqn{n}. If w is not given, all weights are taken to equal 1. Otherwise, the minimization of \eqn{(y - X\beta)'w(y - X\beta)} over \eqn{C} is returned. The default is w = NULL.}
  \item{test}{A logical scalar. If test == TRUE, then the p-value for the test \eqn{H_0:\beta} is in \eqn{V} versus \eqn{H_1:\beta} is in \eqn{C} is returned. \eqn{C} is the constraint cone of the form \eqn{\{\beta: A\beta \ge 0\}}, and \eqn{V} is the null space of \eqn{A}. The default is test = FALSE.}  
  \item{nloop}{ The number of simulations used to get the p-value for the \eqn{E_{01}} test. The default is 1e+4. }
}

\value{
\item{constr.fit }{The constrained fit of \eqn{y} given that \eqn{\beta} is in the cone \eqn{C} of the form \eqn{\{\beta: A\beta \ge 0 \}}.}
\item{unconstr.fit}{The unconstrainted fit, i.e., the least-squares regression of \eqn{y} on the space spanned by \eqn{X}.}
\item{pval}{The p-value for the hypothesis test \eqn{H_0:\beta} is in \eqn{V} versus \eqn{H_1:\beta} is in \eqn{C}. The constraint cone \eqn{C} has the form \eqn{\{\beta: A\beta \ge 0 \}} and \eqn{V} is the null space of \eqn{A}. If test == TRUE, a p-value is returned. Otherwise, the test is skipped and no p-value is returned.}
\item{coefs}{The estimated constrained parameters, i.e., the estimation of the vector \eqn{\beta}.}  
}

\details{
The hypothesis test \eqn{H_0:\beta} is in \eqn{V} versus \eqn{H_1:\beta} is in \eqn{C} is an exact one-sided test, and the test statistic is \eqn{E_{01} = (SSE_0 - SSE_1)/SSE_0}, which has a mixture-of-betas distribution when \eqn{H_0} is true and \eqn{\varepsilon} is a vector following a standard multivariate normal distribution with mean \eqn{0}. The mixing parameters are found through simulations. The number of simulations used to obtain the mixing distribution parameters for the test is 10,000. Such simulations usually take some time. For the "FEV" data set used as an example in this section, whose sample size is 654, the time to get a p-value is roughly 6 seconds.

The constreg function calls coneA for the cone projection part.
}

\references{
Brunk, H. D. (1958) On the estimation of parameters restricted by inequalities.
\emph{The Annals of Mathematical Statistics \bold{29 (2)}}, 437--454.

Raubertas, R. F., C.-I. C. Lee, and E. V. Nordheim (1986) Hypothesis tests for normals
means constrained by linear inequalities. \emph{Communications in Statistics - Theory and 
Methods \bold{15 (9)}}, 2809--2833.

Meyer, M. C. and J. C. Wang (2012) Improved power of one-sided tests. \emph{Statistics and 
Probability Letters \bold{82}}, 1619--1622.  

Liao, X. and M. C. Meyer (2014) coneproj: An R package for the primal or dual cone projections with routines for constrained regression. \emph{Journal of Statistical Software \bold{61(12)}}, 1--22.
}

\author{
Mary C. Meyer and Xiyue Liao 
}

\seealso{
\code{\link{coneA}}
}

\examples{
# load the FEV data set
    data(FEV)

# extract the variables
    y <- FEV$FEV
    age <- FEV$age
    height <- FEV$height
    sex <- FEV$sex
    smoke <- FEV$smoke

# scale age and height
    scale_age <- (age - min(age)) / (max(age) - min(age))
    scale_height <- (height - min(height)) / (max(height) - min(height))

# make xmat
    xmat <- cbind(1, scale_age, scale_height, scale_age * scale_height, sex, smoke)

# make the constraint matrix 
    amat <- matrix(0, 4, 6)
    amat[1, 2] <- 1; amat[2, 2] <- 1; amat[2, 4] <- 1 
    amat[3, 3] <- 1; amat[4, 3] <- 1; amat[4, 4] <- 1

# call constreg to get constrained coefficient estimates
    ans1 <- constreg(y, xmat, amat)
    bhat1 <- coef(ans1)

# call lm to get unconstrained coefficient estimates
    ans2 <- lm(y ~ xmat[,-1])
    bhat2 <- coef(ans2)

# create a 3D plot to show the constrained fit and the unconstrained fit 
    n <- 25
    xgrid <- seq(0, 1, by = 1/n)
    ygrid <- seq(0, 1, by = 1/n)
    x1 <- rep(xgrid, each = length(ygrid))
    x2 <- rep(ygrid, length(xgrid))
    xinterp <- cbind(x1, x2)
    xmatp <- cbind(1, xinterp, x1 * x2, 0, 0)
    
    thint1 <- crossprod(t(xmatp), bhat1)
    A1 <- matrix(thint1, length(xgrid), length(ygrid), byrow = TRUE) 
    thint2 <- crossprod(t(xmatp), bhat2)
    A2 <- matrix(thint2, length(xgrid), length(ygrid), byrow = TRUE) 

    par(mfrow = c(1, 2))
    par(mar = c(4, 1, 1, 1))
    persp(xgrid, ygrid, A1, xlab = "age", ylab = "height", 
    zlab = "FEV", theta = -30)
    title("Constrained Fit")

    par(mar = c(4, 1, 1, 1))
    persp(xgrid, ygrid, A2, xlab = "age", ylab = "height", 
    zlab = "FEV", theta = -30)
    title("Unconstrained Fit")
}
\note{In the 3D plot of the "FEV" example, it is shown that the unconstrained fit increases as "age" increases when "height" is large, but decreases as "age" increases when "height" is small. This does not make sense, since "FEV" should not decrease with respect to "age" given any value of "height". The constrained fit avoids this situation by keeping the fit of "FEV" non-decreasing with respect to "age".
}
\keyword{cone projection routines}
