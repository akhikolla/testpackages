\name{gSOBI}
\alias{gSOBI}
\alias{gSOBI.default}
\alias{gSOBI.ts}
\alias{gSOBI.xts}
\alias{gSOBI.zoo}
\title{
Generalized SOBI
}

\description{
The gSOBI (generalized Second Order Blind Identification) method for the blind source separation (BSS) problem. The method is designed for separating multivariate time series with or without stochastic volatility. The method is a combination of SOBI and vSOBI with \eqn{G(x) = x^2} as a nonlinearity function.
}

\usage{
gSOBI(X, \ldots)

\method{gSOBI}{default}(X, k1 = 1:12, k2 = 1:3, b = 0.9, eps = 1e-06, maxiter = 1000, ordered = FALSE,
      acfk = NULL, original = TRUE, alpha = 0.05, \ldots)
\method{gSOBI}{ts}(X, \ldots)
\method{gSOBI}{xts}(X, \ldots)
\method{gSOBI}{zoo}(X, \ldots)
}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k1}{A vector of lags for SOBI part. It can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:12}.}
  \item{k2}{A vector of lags for vSOBI part. It can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:3}.}
  \item{b}{The weight for the SOBI part, \eqn{1-b} for the vSOBI part. Default is 0.9.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{The maximum number of iterations.}
  \item{ordered}{Whether to order components according to their volatility. Default is \code{FALSE}.}
  \item{acfk}{A vector of lags to be used in testing the presence of serial autocorrelation. Applicable only if \code{ordered = TRUE}.}
  \item{original}{Whether to return the original components or their residuals based on ARMA fit. Default is \code{TRUE}, i.e. the original components are returned. Applicable only if \code{ordered = TRUE}.}
  \item{alpha}{Alpha level for linear correlation detection. Default is 0.05.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
Assume that a \eqn{p}-variate \eqn{{\bf Y}}{Y} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Y}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Y = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{{\bf S}}{S} is the sample covariance matrix of \eqn{{\bf X}}{X}. The algorithm finds an orthogonal matrix \eqn{{\bf U}}{U} by maximizing
\deqn{{\bf D}({\bf U}) = b\sum_{k_1 = 1}^{K_1} {\bf D}_{k_1}({\bf U}) + (1 - b)\sum_{k_2 = 1}^{K_2} {\bf D}_{k_2}({\bf U}),}{D(U) = b*sum_k_1(D_k_1(U)) + (1 - b)*sum_k_2(D_k_2(U)) ,}

where SOBI part

\deqn{{\bf D}_{k_1} = \sum_{i = 1}^p  \left(\frac{1}{T - k_1}\sum_{t=1}^{T - k_1}[({\bf u}_i' {\bf Y}_t) ({\bf u}_i' {\bf Y}_{t + k_1})]\right)^2.}{D_k_1 = sum_i(((1/(T - k_1))*sum_t[(u_i' Y_t)*(u_i' Y_(t + k_1))])^2),}

and vSOBI part

\deqn{{\bf D}_{k_2} = \sum_{i = 1}^p  \left(\frac{1}{T - k_2}\sum_{t=1}^{T - k_2}[({\bf u}_i' {\bf Y}_t)^2 ({\bf u}_i' {\bf Y}_{t + k_2})^2] - \left(\frac{1}{T - k_2}\right)^2\sum_{t=1}^{T - k_2}[({\bf u}_i' {\bf Y}_t)^2]\sum_{t=1}^{T - k_2}[({\bf u}_i' {\bf Y}_{t + k_2})^2]\right)^2}{D_k_2 = sum_i(((1/(T - k_2))*sum_t[(u_i' Y_t)^2*(u_i' Y_(t + k_2))^2] - (1/(T - k_2)^2*sum_t[(u_i' Y_t)^2]*sum_t[(u_i' Y_(t + k_2))^2])^2))}
where \eqn{b \in [0, 1].}{b}
\if{html}{is a value between 0 and 1, and \eqn{i = 1, \ldots, p}, \eqn{k_1 = 1, \ldots, K_1}, \eqn{k_2 = 1, \ldots, K_2} and \eqn{t = 1, \ldots, T}}

The algorithm works iteratively starting with \code{diag(p)} as an initial value for an orthogonal matrix \eqn{{\bf U} = ({\bf u}_1, {\bf u}_2, \ldots, {\bf u}_p)'}{U = (u_1, u_2, \ldots, u_p)'}.

Matrix \eqn{{\bf T}_{ikj}}{T_ikj} is a partial derivative of \eqn{{\bf D}_{kj}({\bf U})}{D_kj(U)}, where \eqn{j = 1, 2}, with respect to \eqn{{\bf u}_i}{u_i}.
Then \eqn{{\bf T}_{kj} = ({\bf T}_{1kj}, \ldots, {\bf T}_{pkj})'}{T_kj = (T_1kj, \ldots, T_pkj)'}, where \eqn{p} is the number of columns in \eqn{\bf Y}{Y}, and \eqn{{\bf T}_j = \sum_{k_j = 1}^{K_j} {\bf T}_{kj}}{T_j = sum(T_kj)}, for \eqn{j = 1, 2}. Finally \eqn{{\bf T} = b{\bf T}_1 + (1-b){\bf T}_2}{T = b*T_1 + (1 - b)*T_2}.

The update for the orthogonal matrix \eqn{{\bf U}_{new} = ({\bf TT}')^{-1/2}{\bf T}}{U.new = (TT')^(-1/2)*T} is calculated at each iteration step. The algorithm stops when
\deqn{||{\bf U}_{new} - {\bf U}_{old}||}{||U.new - U.old||}
is less than \code{eps}.
The final unmixing matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.

For \code{ordered = TRUE} the function orders the sources according to their volatility. First a possible linear autocorrelation is removed using \code{\link[forecast]{auto.arima}}. Then a squared autocorrelation test is performed for the sources (or for their residuals, when linear correlation is present). The sources are then put in a decreasing order according to the value of the test statistic of the squared autocorrelation test. For more information, see \code{\link[tsBSS]{lbtest}}.
}

\value{
 A list with class 'bssvol' containing the following components:
  \item{W }{The estimated unmixing matrix. If \code{ordered = TRUE}, the rows are ordered according to the order of the components.}
  \item{k1}{The vector of the used lags for the SOBI part.}
  \item{k2}{The vector of the used lags for the vSOBI part.}
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

Miettinen, M., Matilainen, M., Nordhausen, K. and Taskinen, S. (2020), \emph{Extracting Conditionally Heteroskedastic Components Using Independent Component Analysis}, Journal of Time Series Analysis,41, 293--311.
}

\author{
Markus Matilainen, Jari Miettinen
}

\seealso{
\code{\link[JADE]{SOBI}}, \code{\link[tsBSS]{vSOBI}}, \code{\link[tsBSS]{lbtest}}, \code{\link[forecast]{auto.arima}}
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

res <- gSOBI(X, 1:4, 1:2, 0.99)
res$W
coef(res)
plot(res)
head(bss.components(res))

MD(res$W, A) # Minimum Distance Index, should be close to zero

# xts series as input
library(xts)
data(sample_matrix)
X2 <- as.xts(sample_matrix)
res2 <- gSOBI(X2, 1:4, 1:2, 0.99)
plot(res2, multi.panel = TRUE)

# zoo series as input
X3 <- as.zoo(X)
res3 <- gSOBI(X3, 1:4, 1:2, 0.99)
plot(res3)
}
}

\keyword{ multivariate }
\keyword{ ts }