\name{tssdr}
\alias{tssdr}
\alias{tssdr.default}
\alias{tssdr.ts}
\alias{tssdr.xts}
\alias{tssdr.zoo}
\alias{print.tssdr}
\alias{components.tssdr}
\alias{plot.tssdr}
\title{
Supervised dimension reduction for multivariate time series
}

\description{
Supervised dimension reduction for multivariate time series data. There are three different algorithms to choose from. TSIR is a time series version of Sliced Inverse Regression (SIR), TSAVE is a time series version of Sliced Average Variance Estimate (TSAVE) and a hybrid of TSIR and TSAVE is TSSH (Time series SIR SAVE Hybrid). For summary of an object of class tssdr, see \code{\link[tsBSS]{summary.tssdr}}.
}

\usage{
tssdr(y, X, \ldots)

\method{tssdr}{default}(y, X, algorithm = c("TSIR", "TSAVE", "TSSH"), k = 1:12, H = 10, weight = 0.5,
      method = c("frjd", "rjd"), eps = 1e-06, maxiter = 1000, \ldots)
\method{tssdr}{ts}(y, X, \ldots)
\method{tssdr}{xts}(y, X, \ldots)
\method{tssdr}{zoo}(y, X, \ldots)

\method{print}{tssdr}(x, digits = 3, ...)
\method{components}{tssdr}(x, ...)
\method{plot}{tssdr}(x, main = "The response and the directions", ...)
}

\arguments{
  \item{y}{A numeric vector or a time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}} (same type as X). Missing values are not allowed.}
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}} (same type as y). Missing values are not allowed.}
  \item{algorithm}{Algorithm to be used. The options are \code{"TSIR"}, \code{"TSAVE"} and \code{"TSSH"}. Default is \code{"TSIR"}.}
  \item{k}{A vector of lags. It can be any non-zero positive integer, or a vector consisting of them. Default is \code{1:12}.}
  \item{H}{The number of slices. If \code{"TSSH"} is used, \eqn{H} is a 2-vector; the first element is used for TSIR part and the second for TSAVE part. Default is \eqn{H = 10}.}
  \item{weight}{Weight \eqn{0 \le a \le 1} for the hybrid method \code{TSSH} only. With \eqn{a = 1} it reduces to TSAVE and with \eqn{a = 0} to TSIR. Default is \eqn{a = 0.5}.}
  \item{method}{The method to use for the joint diagonalization. The options are \code{\link[JADE:rjd]{"rjd"}} and \code{\link[JADE:rjd]{"frjd"}}. Default is \code{\link[JADE:rjd]{"frjd"}}.}
  \item{eps}{Convergence tolerance.}
  \item{maxiter}{The maximum number of iterations.}
  \item{...}{Further arguments to be passed to or from methods.}
  
\strong{In methods for class 'tssdr' only:}

  \item{x}{An object of class tssdr}
  \item{digits}{The number of digits when printing an object of class tssdr. Default is 3}
  \item{main}{A title for a plot when printing an object of class tssdr.}
}

\details{
Assume that the \eqn{p}-variate time series \eqn{{\bf Z}}{Z} with \eqn{T} observations is whitened, i.e. \eqn{{\bf Z}={\bf S}^{-1/2}({\bf X}_t - \frac{1}{T}\sum_{t=1}^T {\bf X}_{t})}{Z = S^(-1/2)*(X_t - (1/T)*sum_t(X_t))}, \if{html}{for \eqn{t = 1, \ldots, T},}
where \eqn{{\bf S}}{S} is a sample covariance matrix of \eqn{{\bf X}}{X}.
Divide \eqn{y} into \eqn{H} disjoint intervals (slices) by its empirical quantiles.

For each lag \eqn{j}, denote \eqn{y_{j}}{y_j} for a vector of the last \eqn{n - j} values of the sliced \eqn{y}. Also denote \eqn{{\bf Z}_j}{Z_j} for the first \eqn{n - j} observations of \eqn{{\bf Z}}{Z}. Then \eqn{{\bf Z}_{jh}}{Z_jh} are the disjoint slices of \eqn{{\bf Z}_j}{Z_j} according to the values of \eqn{y_{j}}{y_j}.

Let \eqn{T_{jh}}{T_jh} be the number of observations in \eqn{{\bf Z}_{jh}}{Z_jh}.
Write \eqn{\bf \widehat{A}_{jh} = \frac{1}{T_{jh}}\sum_{t = 1}^{T_{jh}}({\bf Z}_{jh})_{t}}{A_jh = 1/T_jh * sum_t((Z_jh)_t),} \if{html}{for \eqn{t = 1, \ldots, T_jh},}
and \eqn{{\bf \widehat A}_j = ({\bf \widehat{A}}_{j1}, \ldots, {\bf \widehat{A}}_{jH})'}{Aj = (A_j1, \ldots, A_jH)'}.
Then for algorithm \code{TSIR} matrix \deqn{{\bf \widehat{M}}_{0j} = {\bf \widehat{Cov}}_{A_j}.}{M_0j = Cov_Aj (a sample covariance matrix of Aj).}

Denote \eqn{\bf \widehat{Cov}_{jh}}{Cov_jh} for a sample covariance matrix of \eqn{{\bf Z}_{jh}}{Z_jh}. Then for algorithm \code{TSAVE} matrix \deqn{{\bf \widehat{M}}_{0j} = \frac{1}{H}\sum_{h = 1}^H({\bf I}_p - {\bf \widehat{Cov}_{jh}})^2.}{M_1j = 1/T_jh * sum_h((I_p - Cov_jh)^2),}
\if{html}{\eqn{h = 1, \ldots, H}.}

For \code{TSSH} then matrix \deqn{{\bf \widehat{M}}_{2j} = a{\bf \widehat{M}_{1j}} + (1-a){\bf \widehat{M}_{0j}},}{M_2j = a*M_1j + (1-a)*M_0j,} for a chosen \eqn{0 \le a \le 1}. Note that the value of \eqn{H} can be different for TSIR and TSAVE parts.

The algorithms find an orthogonal matrix \eqn{{\bf U} = (\bf u_1, \ldots, \bf u_p)'}{U = (u_1, \ldots, u_p)'} by maximizing, for \eqn{b = 0, 1} or \eqn{2}, 
\deqn{\sum_{i \in k} ||diag({\bf U} {\bf \widehat{M}}_{bj} {\bf U}')||^2 = \sum_{i \in 1}^p \sum_{j \in k} ({\bf u}_i' {\bf \widehat{M}}_{bj} {\bf u}_i)^2.}{sum_j (||diag(U M_bj U')||^2) = sum_i sum_j (u_i' M_bj u_i)^2,}
\if{html}{for \eqn{i = 1, \ldots, p} and all lags \eqn{j}.}
The final signal separation matrix is then \eqn{{\bf W} = {\bf US}^{-1/2}}{W = U S^(-1/2)}.

Write \eqn{\lambda_{ij} = c({\bf u}_i' {\bf \widehat{M}}_{bj} {\bf u}_i)^2}{\lambda_ij = c*(u_i' M_aj u_i)^2}, where \eqn{c} is chosen in such way that \eqn{\sum_{i = 1}^p \sum_{j \in k} \lambda_{ij}= 1.}{sum_i sum_j \lambda_ij = 1,} \if{html}{for \eqn{i = 1, \ldots, p} and all lags \eqn{j}.}
Then the \eqn{(i, j)}:th element of the matrix \eqn{\bf L}{L} is \eqn{\lambda_{ij}}{\lambda_ij}.

To make a choice on which lags and directions to keep, see \code{\link[tsBSS]{summary.tssdr}}. Note that when printing a tssdr object, all elements are printed, except the directions S.
}

\value{
A list with class 'tssdr' containing the following components:
  \item{W }{The estimated signal separation matrix.}
  \item{k }{The vector of the used lags.}
  \item{S }{The estimated directions as time series object standardized to have mean 0 and unit variances.}
  \item{L }{The Lambda matrix for choosing lags and directions.}
  \item{H }{The used number of slices.}
  \item{yname }{The name for the response time series \eqn{y}.}
  \item{Xname }{The name for the predictor time series \eqn{\bf X}{X}.}
  \item{algorithm }{The used algorithm as a character string.}
}

\references{
Matilainen M., Croux C., Nordhausen K. and Oja H. (2017), \emph{Supervised Dimension Reduction for Multivariate Time Series}, Econometrics and Statistics, 4, 57--69.

Matilainen M., Croux C., Nordhausen K. and Oja H. (2019), \emph{Sliced Average Variance Estimation for Multivariate Time Series}. Statistics: A Journal of Theoretical and Applied Statistics, 53, 630--655, <doi:10.1080/02331888.2019.1605515>.

Li, K.C. (1991), \emph{Sliced Inverse Regression for Dimension Reduction}, Journal of the American Statistical Association, 86, 316--327.

Cook, R. and Weisberg, S. (1991), \emph{Sliced Inverse Regression for Dimension Reduction}, Comment. Journal of the American Statistical Association, 86, 328--332.
}

\author{
Markus Matilainen
}

\seealso{
\code{\link[tsBSS]{summary.tssdr}}, \code{\link[dr]{dr}}
}

\examples{
n <- 10000
A <- matrix(rnorm(9), 3, 3)

x1 <- arima.sim(n = n, list(ar = 0.2))
x2 <- arima.sim(n = n, list(ar = 0.8))
x3 <- arima.sim(n = n, list(ar = 0.3, ma = -0.4))
eps2 <- rnorm(n - 1)
y <- 2*x1[1:(n - 1)] + eps2
X <- ((cbind(x1, x2, x3))[2:n, ]) \%*\% t(A)

res1 <- tssdr(y, X, algorithm = "TSAVE")
res1
summ1 <- summary(res1, type = "alllag", thres = 0.8)
summ1
plot(summ1)
head(components(summ1))
coef(summ1)

# Hybrid of TSIR and TSAVE. For TSIR part H = 10 and for TSAVE part H = 2.
tssdr(y, X, algorithm = "TSSH", weight = 0.6, H = c(10, 2))
}

\keyword{ multivariate }
\keyword{ ts }