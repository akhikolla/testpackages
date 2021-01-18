\name{PVC}
\alias{PVC}
\alias{PVC.default}
\alias{PVC.ts}
\alias{PVC.xts}
\alias{PVC.zoo}
\title{
A modified algorithm for Principal Volatility Component estimator
}

\description{
PVC (Principal Volatility Component) estimator for the blind source separation (BSS) problem. This method is a modified version of PVC by Hu and Tsay (2014).
}

\usage{
PVC(X, \ldots)

\method{PVC}{default}(X, k = 1:12, ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, \ldots)
\method{PVC}{ts}(X, \ldots)
\method{PVC}{xts}(X, \ldots)
\method{PVC}{zoo}(X, \ldots)
}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k}{A vector of lags. It can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:12}.}
  \item{ordered}{Whether to order components according to their volatility. Default is \code{FALSE}.}
  \item{acfk}{A vector of lags to be used in testing the presence of serial autocorrelation. Applicable only if \code{ordered = TRUE}.}
  \item{original}{Whether to return the original components or their residuals based on ARMA fit. Default is \code{TRUE}, i.e. the original components are returned. Applicable only if \code{ordered = TRUE}.}
  \item{alpha}{Alpha level for linear correlation detection. Default is 0.05.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{\bf S}{S} is the sample covariance matrix of \eqn{\bf X}{X}.
Then for each lag \eqn{k} we calculate
\deqn{\widehat{Cov}({\bf Y}_t {\bf Y}_t', Y_{ij, t-k}) = \frac{1}{T}\sum_{t = k + 1}^T \left({\bf Y}_t {\bf Y}_t' - \frac{1}{T-k}\sum_{t = k+1}^T {\bf Y}_t {\bf Y}_t' \right)\left(Y_{ij, t-k} - \frac{1}{T-k}\sum_{t = k+1}^T {Y}_{ij, t-k}\right),}{Cov(Y_t Y_t', Y_(ij, t-k)) = (1/T) * sum_t [{Y_t Y_t' - (1/(T-k))*sum_t [Y_t Y_t']}*{Y_(ij, t - k) - (1/(T-k))*sum_t [Y_(ij, t-k)]}],}
where \if{html}{\eqn{t = k + 1, \ldots, T} and} \eqn{Y_{ij, t-k} = Y_{i, t-k} Y_{j, t-k}, i, j = 1, \ldots, p}{Y_(ij, t-k) = Y_(i, t-k) Y_(j, t-k), i, j = 1, \ldots, p}.
Then 
\deqn{{\bf g}_k({\bf Y}) = \sum_{i = 1}^p \sum_{j=1}^p (\widehat{Cov}({\bf Y}_t {\bf Y}_t', Y_{ij, t-k}))^2.}{g_k(Y) = sum_i [sum_j[{Cov(Y_t Y_t', Y_(ij, t-k))}^2]],}
\if{html}{where \eqn{i,j = 1, \ldots, p.}}
Thus the generalized kurtosis matrix is 
\deqn{{\bf G}_K({\bf Y}) = \sum_{k = 1}^K {\bf g}_k({\bf Y}),}{G_K(Y) = sum_k [g_k(Y)],}
where \eqn{k = 1, \ldots, K} is the set of chosen lags.
Then \eqn{\bf U}{U} is the matrix with eigenvectors of \eqn{{\bf G}_K({\bf Y})}{G_K(Y)} as its rows.
The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}, where the average value of each row is set to be positive.

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
Miettinen, M., Matilainen, M., Nordhausen, K. and Taskinen, S. (2020), \emph{Extracting Conditionally Heteroskedastic Components Using Independent Component Analysis}, Journal of Time Series Analysis,41, 293--311.

Hu, Y.-P. and Tsay, R. S. (2014), \emph{Principal Volatility Component Analysis}, Journal of Business & Economic Statistics, 32(2), 153--164.
}

\author{
Jari Miettinen, Markus Matilainen
}

\seealso{
\code{\link[MTS]{comVol}}, \code{\link[tsBSS]{gSOBI}}, \code{\link[tsBSS]{lbtest}}, \code{\link[forecast]{auto.arima}}
}

\examples{
if(require(stochvol)) {
n <- 10000
A <- matrix(rnorm(9), 3, 3)

# simulate SV models
s1 <- svsim(n, mu = -10, phi = 0.8, sigma = 0.1)$y
s2 <- svsim(n, mu = -10, phi = 0.9, sigma = 0.2)$y
s3 <- svsim(n, mu = -10, phi = 0.95, sigma = 0.4)$y

# create a daily time series
X <- ts(cbind(s1, s2, s3) \%*\% t(A), end = c(2015, 338), frequency = 365.25)

res <- PVC(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero
}
}

\keyword{ multivariate }
\keyword{ ts }