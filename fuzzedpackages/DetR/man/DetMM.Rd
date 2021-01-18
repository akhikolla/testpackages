\name{DetMM}
\alias{DetMM}
\title{Robust and Deterministic Linear Regression via DetMM}
\description{
     Function to compute the DetMM estimates of regression.
}
\usage{
  DetMM(x,y,intercept=1,alpha=0.75,h=NULL,scale_est="scaleTau2",tuning.chi=1.54764,
tuning.psi=4.685061)
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
\item{tuning.chi}{tuning constant vector for the bi-weight chi used for the ISteps.}
\item{tuning.psi}{tuning constant vector for the bi-weight psi used for the MSteps.}
}
\value{
  The function \code{DetLTS} returns a list with as many components as 
	there are elements in the \code{h}. Each of the entries is a list
containing  the following components:
	\item{coefficients}{The estimate of the coefficient vector}
  \item{scale}{The scale as used in the M steps.}
  \item{residuals}{Residuals associated with the estimator.}
  %loss
  \item{converged}{\code{TRUE} if the IRWLS iterations have converged.}
  \item{iter}{number of IRWLS iterations}
  \item{rweights}{the \dQuote{robustness weights} \eqn{\psi(r_i/S) / (r_i/S)}.}
  \item{fitted.values}{Fitted values associated with the estimator.}
  \item{DetS}{A similar list that contains the results of (initial) returned by DetS}
}
\references{
    Maronna, R.A. and Zamar, R.H. (2002) Robust estimates of location
     and dispersion of high-dimensional datasets; \emph{Technometrics}
     \bold{44}(4), 307--317.

 Rousseeuw, P.J. and Croux, C. (1993) Alternatives to the Median Absolute Deviation; 
  \emph{Journal of the American Statistical Association} , \bold{88}(424), 1273--1283.

  Croux, C., Dhaene, G. and Hoorelbeke, D. (2003)
  \emph{Robust standard errors for robust estimators},
  Discussion Papers Series 03.16, K.U. Leuven, CES.

  Koller, M. (2012), Nonsingular subsampling for S-estimators with
  categorical predictors, \emph{ArXiv e-prints}, arXiv:1208.5595v1.

  Koller, M. and Stahel, W.A. (2011), Sharpening Wald-type inference in
  robust regression for small samples, \emph{Computational Statistics &
    Data Analysis} \bold{55}(8), 2504--2515.

  Maronna, R. A., and Yohai, V. J. (2000).
  Robust regression with both continuous and categorical predictors.
  \emph{Journal of Statistical Planning and Inference} \bold{89}, 197--214.

  Rousseeuw, P.J. and Yohai, V.J. (1984)
  Robust regression by means of S-estimators,
  In \emph{Robust and Nonlinear Time Series},
  J. Franke, W. Hardle and R. D. Martin (eds.).
  Lectures Notes in Statistics 26, 256--272,
  Springer Verlag, New York.

  Salibian-Barrera, M. and Yohai, V.J. (2006)
  A fast algorithm for S-regression estimates,
  \emph{Journal of Computational and Graphical Statistics},
  \bold{15}(2), 414--427.

  Yohai, V.J. (1987)
  High breakdown-point and high efficiency estimates for regression.
  \emph{The Annals of Statistics} \bold{15}, 642--65.

}
\examples{
## generate data
set.seed(1234)  # for reproducibility
n<-100
h<-c(55,76,89)
set.seed(123)
x0<-matrix(rnorm(n*2),nc=2)
y0<-rnorm(n)
out1<-DetMM(x0,y0,h=h)
}
\author{
 Vakili Kaveh using translation of the C code from pcaPP (by Peter Filzmoser, Heinrich Fritz, Klaudius Kalcher, see citation("pcaPP")) for 
 the Qn and scaleTau2 (Original by Kjell Konis with substantial modifications by Martin Maechler)  from robustbase
 (see citation("scaleTau2")). This function calls lmrob in package robustbase.
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
