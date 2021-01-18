\name{lbtest}
\alias{lbtest}
\alias{print.lbtest}
\title{
Modified Ljung-Box test and volatility clustering test for time series.
}

\description{
Modified Ljung-Box test and volatility clustering test for time series. Time series can be univariate or multivariate. The modified Ljung-Box test checks whether there is linear autocorrelation in the time series. The volatility clustering test checks whether the time series has squared autocorrelation, which would indicate a presence of volatility clustering.
}

\usage{
lbtest(X, k, type = c("squared", "linear"))

\method{print}{lbtest}(x, digits = 3, ...)
}

\arguments{
  \item{X}{A numeric vector/matrix or a univariate/multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k}{A vector of lags.}
  \item{type}{The type of the autocorrelation test. Options are Modified Ljung-Box test (\code{"linear"}) or volatility clustering test (\code{"squared"}) autocorrelation. Default is \code{"squared"}.}
In methods for class 'lbtest' only:
  \item{x}{An object of class lbtest}
  \item{digits}{The number of digits when printing an object of class lbtest. Default is 3}
  \item{\dots}{Further arguments to be passed to or from methods.}
}

\details{
Assume all the individual time series \eqn{X_i} in \eqn{\bf X}{X} with \eqn{T} observations are scaled to have variance 1.

Then the modified Ljung-Box test statistic for testing the existence of linear autocorrelation in \eqn{X_i} (\code{option = "linear"}) is
\deqn{T \sum_{j \in k} \left(\sum_{t=1}^T (X_{it} X_{i, t + j})/(T - j)\right)^2/V_{j}.}{T*sum_j(sum_t (X_it X_(i, t + j))/(T - j))^2/V_j.}

Here  
\deqn{V_{j} = \sum_{t=1}^{n-j}\frac{x_t^2 x_{t+j}^2}{n-j} + 2 \sum_{k=1}^{n-j-1} \frac{n-k}{n} \sum_{s=1}^{n-k-j}\frac{x_s x_{s+j }x_{s+k} x_{s+k+j}}{n-k-j}.}{V_j = sum_t (x_t^2 x_(t+j)^2)/(n-j) + 2*sum_k (n - k)/n * sum_s (x_s x_(s+j) x_(s+k) x_(s+k+j))/(n - k - j),}
\if{html}{where \eqn{t = 1, \ldots, n - j}, \eqn{k = 1, \ldots, n - j - 1} and \eqn{s = 1, \ldots, n - k - j}.}

The volatility clustering test statistic (\code{option = "squared"}) is
\deqn{T \sum_{j \in k} \left(\sum_{t=1}^T (X_{it}^2 X_{i, t + j}^2)/(T - j) - 1\right)^2}{T*sum_j(sum_t (X_it^2 X_(i, t + j)^2)/(T - j) - 1)^2/4.}

Test statistic related to each time series \eqn{X_i} is then compared to \eqn{\chi^2}-distribution with \code{length(k)} degrees of freedom, and the corresponding p-values are produced. Small p-value indicates the existence of autocorrelation.
}

\value{
 A list with class 'lbtest' containing the following components:
  \item{TS}{The values of the test statistic for each component of X as a vector.}
  \item{p_val}{The p-values based on the test statistic for each component of X as a vector.}
  \item{Xname}{The name of the data used as a character string.}
  \item{varnames}{The names of the variables used as a character string vector.}
  \item{k}{The lags used for testing the serial autocorrelation as a vector.}
  \item{K}{The total number of lags used for testing the serial autocorrelation.}
  \item{type}{The type of the autocorrelation test.}
}

\references{
Miettinen, M., Matilainen, M., Nordhausen, K. and Taskinen, S. (2020), \emph{Extracting Conditionally Heteroskedastic Components Using Independent Component Analysis}, Journal of Time Series Analysis,41, 293--311.
}

\author{
Markus Matilainen, Jari Miettinen
}

\seealso{
\code{\link[tsBSS]{FixNA}}, \code{\link[tsBSS]{gFOBI}}, \code{\link[tsBSS]{gJADE}}, \code{\link[tsBSS]{vSOBI}}, \code{\link[tsBSS]{gSOBI}}
}

\examples{
if(require(stochvol)) {
n <- 10000
s1 <- svsim(n, mu = -10, phi = 0.95, sigma = 0.1)$y
s2 <- rnorm(n)
S <- cbind(s1, s2)

lbtest(S, 1:3, type = "squared")
# First p-value should be very close to zero, as there exists stochastic volatility
}
}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ ts }