\name{vSOBI}
\alias{vSOBI}
\alias{vSOBI.default}
\alias{vSOBI.ts}
\alias{vSOBI.xts}
\alias{vSOBI.zoo}
\title{
A Variant of SOBI for Blind Source Separation
}

\description{
The vSOBI (variant of Second Order Blind Identification) method for the blind source separation of time series with stochastic volatility.
The method is a variant of SOBI, which is a method designed to separate ARMA sources, and an alternative to FixNA and FixNA2 methods.
}

\usage{
vSOBI(X, \ldots)

\method{vSOBI}{default}(X, k = 1:12, eps = 1e-06, maxiter = 1000, G = c("pow", "lcosh"), 
                        ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, \ldots)
\method{vSOBI}{ts}(X, \ldots)
\method{vSOBI}{xts}(X, \ldots)
\method{vSOBI}{zoo}(X, \ldots)
}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k}{A vector of lags. It can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:12}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{The maximum number of iterations.}
  \item{G}{Function \eqn{G(x)}. The choices are \code{"pow"} (default) and \code{"lcosh"}.}
  \item{ordered}{Whether to order components according to their volatility. Default is \code{FALSE}.}
  \item{acfk}{A vector of lags to be used in testing the presence of serial autocorrelation. Applicable only if \code{ordered = TRUE}.}
  \item{original}{Whether to return the original components or their residuals based on ARMA fit. Default is \code{TRUE}, i.e. the original components are returned. Applicable only if \code{ordered = TRUE}.}
  \item{alpha}{Alpha level for linear correlation detection. Default is 0.05.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{{\bf S}}{S} is the sample covariance matrix of \eqn{{\bf X}}{X}. The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing
\deqn{{\bf D}({\bf U}) = \sum_{k = 1}^K {\bf D}_k({\bf U})}{D(U) = sum_k(D_k(U))}
\deqn{= \sum_{k = 1}^K \sum_{i = 1}^p  \left(\frac{1}{T - k}\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_t) G({\bf u}_i' {\bf Y}_{t + k})] - \left(\frac{1}{T - k}\right)^2\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_t)]\sum_{t=1}^{T - k}[G({\bf u}_i' {\bf Y}_{t + k})]\right)^2.}{ = sum_k(sum_i(((1/(T - k))*sum_t[G(u_i' Y_t)*G(u_i' Y_(t + k))] - (1/(T - k)^2*sum_t[G(u_i' Y_t)]*sum_t[G(u_i' Y_(t + k))])^2))),}
\if{html}{where \eqn{i = 1, \ldots, p}, \eqn{k = 1, \ldots, K} and \eqn{t = 1, \ldots, T}.} For function \eqn{G(x)} the choices are \eqn{x^2} and log(cosh(\eqn{x})).

The algorithm works iteratively starting with \code{diag(p)} as an initial value for an orthogonal matrix \eqn{{\bf U} = ({\bf u}_1, {\bf u}_2, \ldots, {\bf u}_p)'}{U = (u_1, u_2, \ldots, u_p)'}.

Matrix \eqn{{\bf T}_{ik}}{T_ik} is a partial derivative of \eqn{{\bf D}_k({\bf U})}{D_k(U)} with respect to \eqn{{\bf u}_i}{u_i}.
Then \eqn{{\bf T}_k = ({\bf T}_{1k}, \ldots, {\bf T}_{pk})'}{T_k = (T_1k, \ldots, T_pk)'}, where \eqn{p} is the number of columns in \eqn{{\bf Y}}{Y}, and \eqn{{\bf T} = \sum_{k = 1}^K {\bf T}_k}{T = sum(T_k)}.
The update for the orthogonal matrix \eqn{{\bf U}_{new} = ({\bf TT}')^{-1/2}{\bf T}}{U.new = (TT')^(-1/2)*T} is calculated at each iteration step. The algorithm stops when
\deqn{||{\bf U}_{new} - {\bf U}_{old}||}{||U.new - U.old||}
is less than \code{eps}.
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
Belouchrani, A., Abed-Meriam, K., Cardoso, J.F. and Moulines, R. (1997), \emph{A Blind Source Separation Technique Using Second-Order Statistics}, IEEE Transactions on Signal Processing, 434--444.

Matilainen, M., Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2017), \emph{On Independent Component Analysis with Stochastic Volatility Models}, Austrian Journal of Statistics, 46(3--4), 57--66.
}

\author{
Markus Matilainen
}

\seealso{
\code{\link[tsBSS]{FixNA}}, \code{\link[JADE]{SOBI}}, \code{\link[tsBSS]{lbtest}}, \code{\link[forecast]{auto.arima}}
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

res <- vSOBI(X)
res
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero
}
}

\keyword{ multivariate }
\keyword{ ts }