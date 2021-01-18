\name{gFOBI}
\alias{gFOBI}
\alias{gFOBI.default}
\alias{gFOBI.ts}
\alias{gFOBI.xts}
\alias{gFOBI.zoo}
\title{
Generalized FOBI
}

\description{
The gFOBI (generalized Fourth Order Blind Identification) method for blind source separation of time series with stochastic volatility. The method is a generalization of FOBI, which is a method designed for iid data.

}

\usage{
gFOBI(X, \ldots)

\method{gFOBI}{default}(X, k = 0:12, eps = 1e-06, maxiter = 100, method = c("frjd", "rjd"),
                        na.action = na.fail, weight = NULL, ordered = FALSE,
                        acfk = NULL, original = TRUE, alpha = 0.05, \ldots)
\method{gFOBI}{ts}(X, \ldots)
\method{gFOBI}{xts}(X, \ldots)
\method{gFOBI}{zoo}(X, \ldots)
}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k}{A vector of lags. It can be any non-negative integer, or a vector consisting of them. Default is \code{0:12}. If \eqn{k = 0}, this method reduces to \code{\link[JADE]{FOBI}}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{The maximum number of iterations.}
  \item{method}{The method to use for the joint diagonalization. The options are \code{\link[JADE:rjd]{"rjd"}} and \code{\link[JADE:rjd]{"frjd"}}. Default is \code{\link[JADE:rjd]{"frjd"}}.}
  \item{na.action}{A function which indicates what should happen when the data contain 'NA's. Default is to fail.}
  \item{weight}{A vector of length k to give weight to the different matrices in joint diagonalization. If NULL, all matrices have equal weight.}
  \item{ordered}{Whether to order components according to their volatility. Default is \code{FALSE}.}
  \item{acfk}{A vector of lags to be used in testing the presence of serial autocorrelation. Applicable only if \code{ordered = TRUE}.}
  \item{original}{Whether to return the original components or their residuals based on ARMA fit. Default is \code{TRUE}, i.e. the original components are returned. Applicable only if \code{ordered = TRUE}.}
  \item{alpha}{Alpha level for linear correlation detection. Default is 0.05.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{{\bf S}}{S} is the sample covariance matrix of \eqn{{\bf X}}{X}. Algorithm first calculates
\deqn{{\bf \widehat{B}}^{ij}_k({\bf Y}) =  \frac{1}{T - k} \sum_{t = 1}^T [{\bf Y}_{t + k} {\bf Y}_t' {\bf E}^{ij} {\bf Y}_t {\bf Y}_{t + k}']}{B^ij_k(Y) = (1/(T - k))*sum[Y_(t + k) Y_t'  E^ij Y_t Y_(t + k)'],}
\if{html}{where \eqn{t = 1, \ldots, T},} and then 
\deqn{{\bf \widehat{B}}_k({\bf Y}) = \sum_{i = 1}^p {\bf \widehat{B}}^{ii}_k({\bf Y}).}{B_k(Y) = sum(B^ii_k(Y)),}
\if{html}{for \eqn{i = 1, \ldots, p}.}

The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing 
\deqn{\sum_{k = 0}^K ||\textrm{diag}({\bf U \widehat{B}}_k({\bf Y}) {\bf U}')||^2.}{sum(||diag(U B_k(Y) U')||^2).}
The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.

For \code{ordered = TRUE} the function orders the sources according to their volatility. First a possible linear autocorrelation is removed using \code{\link[forecast]{auto.arima}}. Then a squared autocorrelation test is performed for the sources (or for their residuals, when linear correlation is present). The sources are then put in a decreasing order according to the value of the test statistic of the squared autocorrelation test. For more information, see \code{\link[tsBSS]{lbtest}}.
}

\value{
 A list with class 'bssvol' containing the following components:
  \item{W }{The estimated unmixing matrix. If \code{ordered = TRUE}, the rows are ordered according to the order of the components.}
  \item{k }{The vector of the used lags.}
  \item{S }{The estimated sources as time series object standardized to have mean 0 and unit variances. If \code{ordered = TRUE}, then components are ordered according to their volatility. If \code{original = FALSE}, the sources with linear autocorrelation are replaced by their ARMA residuals.}
  \item{MU}{The means of the original series.}
  If \code{ordered = TRUE}, then also the following components included in the list:
  \item{Sraw }{The ordered original estimated sources as time series object standardized to have mean 0 and unit variances. Returned only if \code{original = FALSE}.}
  \item{fits }{The ARMA fits for the components with linear autocorrelation.}
  \item{armaeff }{A logical vector. Is TRUE if ARMA fit was done to the corresponding component.}
  \item{linTS }{The value of the modified Ljung-Box test statistic for each component.}
  \item{linP }{p-value based on the modified Ljung-Box test statistic for each component.}
  \item{volTS }{The value of the volatility clustering test statistic.}
  \item{volP }{p-value based on the volatility clustering test statistic.}
}

\references{
Cardoso, J.-F. (1989), \emph{Source Separation Using Higher Order Moments}, in: Proceedings of the IEEE International Conference on Acoustics, Speech and Signal Processing, 2109--2112.

Matilainen, M., Nordhausen, K. and Oja, H. (2015), \emph{New Independent Component Analysis Tools for Time Series}, Statistics & Probability Letters, 105, 80--87.
}

\author{
Markus Matilainen, Klaus Nordhausen
}

\seealso{
\code{\link[JADE]{FOBI}}, \code{\link[JADE:rjd]{frjd}}, \code{\link[tsBSS]{lbtest}}, \code{\link[forecast]{auto.arima}}
}

\examples{
if(require(stochvol)) {
n <- 10000
A <- matrix(rnorm(9), 3, 3)

# simulate SV models
s1 <- svsim(n, mu = -10, phi = 0.8, sigma = 0.1)$y
s2 <- svsim(n, mu = -10, phi = 0.9, sigma = 0.2)$y
s3 <- svsim(n, mu = -10, phi = 0.95, sigma = 0.4)$y

X <- cbind(s1, s2, s3) \%*\% t(A)

res <- gFOBI(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero
}
}

\keyword{ multivariate }
\keyword{ ts }