\name{DetLTS}
\alias{DetLTS}
\title{Robust and Deterministic Linear Regression via DetLTS}
\description{
     Function to compute the DetLTS estimates of regression.
}
\usage{
  DetLTS(x, y, intercept = 1, alpha = 0.75, h = NULL, scale_est = "scaleTau2")
}
\arguments{
  \item{x}{Matrix of design variables. Never contains an intercept.}
  \item{y}{Vector of responses.}
  \item{intercept}{A boolean indicating whether the regression contains an intercept.}
  \item{alpha}{numeric parameter controlling the size of the subsets over which the determinant is minimized, i.e., alpha*n observations are used for computing the determinant. Allowed values are between 0.5 and 1 and the default is 0.75. Can be a vector.}
\item{h}{Integer in [\code{ceiling((n+p+1)/2)},n) which determines the number of observations which 
are awarded weight in the fitting process. Can be a vector. If both \code{h} and \code{alpha} are set to non default values,
 \code{alpha} will be ignored.}
  \item{scale_est}{A character string specifying the
  variance functional.  Possible values are "Qn" or "scaleTau2".}
}
\value{
  The function \code{DetLTS} returns a list with as many components as 
	there are elements in the \code{h}. Each of the entries is a list
containing  the following components:
  \item{crit}{
    the value of the objective function of the LTS regression method,
    i.e., the sum of the \eqn{h} smallest squared raw residuals.
  }
  \item{coefficients}{
    vector of coefficient estimates (including the intercept by default when
    \code{intercept=TRUE}), obtained after reweighting.
  }
  \item{best}{
    the best subset found and used for computing the raw estimates, with
    \code{\link{length}(best) == quan = \link{h.alpha.n}(alpha,n,p)}.
  }
  \item{fitted.values}{vector like \code{y} containing the fitted values
    of the response after reweighting.}
  \item{residuals}{vector like \code{y} containing the residuals from
    the weighted least squares regression.}
  \item{scale}{scale estimate of the reweighted residuals.  }
  \item{alpha}{same as the input parameter \code{alpha}.}
  \item{quan}{the number \eqn{h} of observations which have determined
    the least trimmed squares estimator.}
  \item{intercept}{same as the input parameter \code{intercept}.}
  \item{cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of
    the final estimate of the error scale.}
  \item{raw.coefficients}{vector of raw coefficient estimates (including
  the intercept, when \code{intercept=TRUE}).}
  \item{raw.scale}{scale estimate of the raw residuals.}
  \item{raw.resid}{vector like \code{y} containing the raw residuals
    from the regression.}
  \item{raw.cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of the
    raw estimate of the error scale.}
  \item{lts.wt}{
    vector like y containing weights that can be used in a weighted
    least squares.  These weights are 1 for points with reasonably
    small residuals, and 0 for points with large residuals.
  }
  \item{raw.weights}{
    vector containing the raw weights based on the raw residuals and raw scale.
  }
  \item{method}{character string naming the method (Least Trimmed Squares).}

}
\references{

    Vakili K. (2016). A study and implementation of robust estimators for 
    multivariate and functional data (Doctoral dissertation). 	

    Maronna, R.A. and Zamar, R.H. (2002) Robust estimates of location
     and dispersion of high-dimensional datasets; \emph{Technometrics}
     \bold{44}(4), 307--317.

 Rousseeuw, P.J. and Croux, C. (1993) Alternatives to the Median Absolute Deviation; 
  \emph{Journal of the American Statistical Association} , \bold{88}(424), 1273--1283.

  Peter J. Rousseeuw (1984), Least Median of Squares Regression.
  \emph{Journal of the American Statistical Association} \bold{79}, 871--881.

  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.

  Pison, G., Van Aelst, S., and Willems, G. (2002)
  Small Sample Corrections for LTS and MCD.
  \emph{Metrika} \bold{55}, 111-123.
}
\examples{
n<-100
h<-c(55,76,89)
set.seed(123)# for reproducibility
x0<-matrix(rnorm(n*2),nc=2)
y0<-rnorm(n)
out1<-DetLTS(x0,y0,h=h)
}
\author{
 Vakili Kaveh using translation of the C code from pcaPP (by Peter Filzmoser, Heinrich Fritz, Klaudius Kalcher, see citation("pcaPP")) for 
 the Qn and scaleTau2 (Original by Kjell Konis with substantial modifications by Martin Maechler)  from robustbase
 (see citation("scaleTau2"))
 as well as R code from function ltsReg in package robustbase (originally written by Valentin Todorov valentin.todorov@chello.at, based on work written for S-plus by Peter Rousseeuw and Katrien van Driessen from
  University of Antwerp, see citation("ltsReg")).
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
