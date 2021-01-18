\name{summary.tssdr}
\alias{summary.tssdr}
\alias{print.summary.tssdr}
\alias{components.summary.tssdr}
\alias{coef.summary.tssdr}
\alias{plot.summary.tssdr}
\title{
Summary of an object of class tssdr
}

\description{
Gives a summary of an object of class tssdr. It includes different types of methods to select the number of directions (sources) and lags.
}

\usage{
\method{summary}{tssdr}(object, type = c("rectangle", "alllag", "alldir", "big"), thres = 0.8, ...)

\method{print}{summary.tssdr}(x, digits = 3, ...)
\method{components}{summary.tssdr}(x, ...)
\method{coef}{summary.tssdr}(object, ...)
\method{plot}{summary.tssdr}(x, main = "The response and the chosen directions", ...)
}

\arguments{
  \item{object}{An object of class tssdr.}
  \item{type}{Method for choosing the important lags and directions. The choices are \code{"rectangle"}, \code{"alllag"}, \code{"alldir"} and \code{"big"}. Default is \code{"rectangle"}.}
  \item{thres}{The threshold value for choosing the lags and directions. Default is \code{0.8}.}
  \item{\dots}{Further arguments to be passed to or from methods.}
  
\strong{In methods for class 'summary.tssdr' only:}
  \item{x}{An object of class summary.tssdr}
  \item{digits}{The number of digits when printing an object of class summary.tssdr. Default is 3}
  \item{main}{A title for a plot when printing an object of class summary.tssdr.}
}

\details{
The sum of values of \eqn{k_0 \times p_0}{k0*p0} matrix \eqn{\bf L}{L} of \code{object} is 1. The values of the matrix are summed together in ways detailed below, until the value is at least \eqn{\pi} (\code{thres}). Let \eqn{\lambda_{ij}}{\lambda_ij} be the element \eqn{(i, j)} of the matrix \eqn{\bf L}{L}.

For \code{alllag}: \eqn{k = k_0}{k = k0} and \eqn{p} is the smallest value for which \eqn{\sum_{i = 1}^p \lambda_{ij} \ge \pi.}{sum_i (\lambda_ij) \ge \pi,} \if{html}{where \eqn{i = 1, \ldots, p}. The chosen number of lags and directions are returned.}

For \code{alldir}: \eqn{p = p_0}{p = p0} and \eqn{k} is the smallest value for which \eqn{\sum_{j = 1}^k \lambda_{ij} \ge \pi}{sum_j (\lambda_ij) \ge \pi,} \if{html}{where \eqn{j = 1, \ldots, k}. The chosen number of lags and directions are returned.}

For \code{rectangle}: \eqn{k} and \eqn{p} are values such that their product \eqn{k p}{ k*p} is the smallest for which \eqn{\sum_{i = 1}^p \sum_{j = 1}^k \lambda_{ij} \ge \pi}{sum_i sum_j (\lambda_ij) \ge \pi,} \if{html}{where \eqn{i = 1, \ldots, p} and \eqn{j = 1, \ldots, k}. The chosen number of lags and directions are returned.}


For \code{big}: \eqn{r} is the smallest value of elements \eqn{(i_1, j_1), \ldots, (i_r, j_r)} for which \eqn{\sum_{k = 1}^r \lambda_{i_k,j_k} \ge \pi}{sum_k(\lambda_(i_r,j_r)) \ge \pi,} \if{html}{where \eqn{k = 1, \ldots, r}. Thi indices of the matrix corresponding to the chosen values are returned.}

Note that when printing a summary.tssdr object, all elements except the component S, which is the matrix of the chosen directions or a vector if there is only one direction, are printed.
}

\value{
A list with class 'summary.tssdr' containing the following components:
  \item{W }{The estimated signal separation matrix}
  \item{L }{The Lambda matrix for choosing lags and directions.}
  \item{S }{The estimated directions as time series object standardized to have mean 0 and unit variances.}
  \item{type}{The method for choosing the important lags and directions.}
  \item{algorithm }{The used algorithm as a character string.}
  \item{yname }{The name for the response time series \eqn{y}.}
  \item{Xname }{The name for the predictor time series \eqn{\bf X}{X}.}
  \item{k }{The chosen number of lags (not for \code{type = "big"} ).}
  \item{p }{The chosen number of directions (not for \code{type = "big"}).}
  \item{pk}{The chosen lag-direction combinations (for \code{type = "big"} only).}
}

\references{
Matilainen M., Croux C., Nordhausen K. and Oja H. (2017), \emph{Supervised Dimension Reduction for Multivariate Time Series}, Econometrics and Statistics, 4, 57--69.
}

\author{Markus Matilainen}

\seealso{
\code{\link[tsBSS]{tssdr}}
}

\examples{
n <- 10000
A <- matrix(rnorm(9), 3, 3)

x1 <- arima.sim(n = n, list(ar = 0.2))
x2 <- arima.sim(n = n, list(ar = 0.8))
x3 <- arima.sim(n = n, list(ar = 0.3, ma = -0.4))
eps2 <- rnorm(n - 1)
y <- 2*x1[1:(n - 1)] + 3*x2[1:(n - 1)] + eps2
X <- ((cbind(x1, x2, x3))[2:n, ]) \%*\% t(A)

res2 <- tssdr(y, X, algorithm = "TSIR")
res2
summ2 <- summary(res2, thres = 0.5)
summ2
summary(res2) #Chooses more lags with larger threshold
summary(res2, type = "alllag") #Chooses all lags
summary(res2, type = "alldir", thres = 0.5) #Chooses all directions
summary(res2, type = "big", thres = 0.5) #Same choices than in summ2
}

\keyword{methods}