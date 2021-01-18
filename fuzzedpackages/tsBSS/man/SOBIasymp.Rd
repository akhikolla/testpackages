\name{SOBIasymp}
\alias{SOBIasymp}
\alias{SOBIasymp.default}
\alias{SOBIasymp.ts}
\alias{SOBIasymp.xts}
\alias{SOBIasymp.zoo}
\title{
Second-order Separation Sub-White-Noise Asymptotic Testing with SOBI
}

\description{
The function uses SOBI (Second Order Blind Identification) to test whether the last \code{p-k} latent series are pure white noise, assuming a p-variate second-order stationary blind source separation (BSS) model. The test is asymptotic.
}

\usage{
SOBIasymp(X, \ldots)

\method{SOBIasymp}{default}(X, k, tau = 1:12, eps = 1e-06, maxiter = 200, \ldots)
\method{SOBIasymp}{ts}(X, \ldots)
\method{SOBIasymp}{xts}(X, \ldots)
\method{SOBIasymp}{zoo}(X, \ldots)

}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k}{The number of latent series that are not white noise. Can be between \eqn{0} and \eqn{p-1}.}
  \item{tau}{The lags for the SOBI autocovariance matrices.}
  \item{eps}{The convergence tolerance for the joint diagonalization.}
  \item{maxiter}{The maximum number of iterations for the joint diagonalization.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
SOBI standardizes \code{X} with \eqn{n} samples and jointly diagonalizes the autocovariance matrices of the standardized data for a chosen set of lags \code{tau}, yielding a transformation \eqn{\bf W}{W} giving the latent variables as \eqn{{\bf S} = {\bf X} {\bf W}}{S = XW}. Assume, without loss of generality, that the latent components are ordered in decreasing order with respect to the sums of squares of the corresponding "eigenvalues" produced by the joint diagonalization. Under the null hypothesis the lower right corner \eqn{(p - k) \times (p - k)}{(p - k)*(p - k)} blocks of the autocovariance matrices of the sources are zero matrices and the sum \eqn{m} of their squared norms over all lags can be used as a test statistic in inference on the true number of latent white noise series.

This function conducts the hypothesis test using the asymptotic null distribution of \eqn{m}, a chi-squared distribution with \eqn{T(p - k)(p - k + 1)/2} degrees of freedom where \eqn{T} is the number of autocovariance matrices used by SOBI.
}

\value{
A list of class ictest inheriting from class htest containing:
  \item{statistic}{The value of the test statistic.}
  \item{p.value}{The p-value of the test.}
  \item{parameter}{The degrees of freedom of the asymptotic null distribution.}
  \item{method}{Character string indicating which test was performed.}
  \item{data.name}{Character string giving the name of the data.}
  \item{alternative}{Character string specifying the alternative hypothesis.}
  \item{k}{The number of latent series that are not white noise used in the testing problem.}
  \item{W}{The transformation matrix to the latent series.}
  \item{S}{Multivariate time series with the centered source components.}
  \item{D}{The underlying eigenvalues of the autocovariance matrix.}
  \item{MU}{The location of the data which was subtracted before calculating SOBI.}
  \item{tau}{The used set of lags for the SOBI autocovariance matrices.}
}

\references{
Virta, J. and Nordhausen, K. (2019), \emph{Determining the Signal Dimension in Second Order Source Separation}. To appear in Statistica Sinica, <doi:10.5705/ss.202018.0347>.
}

\author{
Klaus Nordhausen, Joni Virta
}

\seealso{
\code{\link[JADE]{AMUSE}}, \code{\link[JADE]{SOBI}}, \code{\link{AMUSEasymp}}
}

\examples{
  n <- 1000

  A <- matrix(rnorm(16), 4, 4)
  s1 <- arima.sim(list(ar = c(0, 0.6)), n)
  s2 <- arima.sim(list(ma = c(0, -0.5)), n)
  s3 <- rnorm(n)
  s4 <- rnorm(n)

  S <- cbind(s1, s2, s3, s4)
  X <- S \%*\% t(A)

  asymp_res_1 <- SOBIasymp(X, k = 1)
  asymp_res_1

  asymp_res_2 <- SOBIasymp(X, k = 2)
  asymp_res_2

  # Plots of the estimated sources, the last two are white noise
  plot(asymp_res_2)
  
  # Note that AMUSEasymp with lag 1 does not work due to the lack of short range dependencies
  AMUSEasymp(X, k = 1)
}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ ts }