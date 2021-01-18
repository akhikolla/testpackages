\name{AMUSEasymp}
\alias{AMUSEasymp}
\alias{AMUSEasymp.default}
\alias{AMUSEasymp.ts}
\alias{AMUSEasymp.xts}
\alias{AMUSEasymp.zoo}
\title{
Second-order Separation Sub-White-Noise Asymptotic Testing with AMUSE
}

\description{
The function uses AMUSE (Algorithm for Multiple Unknown Signals Extraction) to test whether the last \code{p-k} latent series are pure white noise, assuming a p-variate second-order stationary blind source separation (BSS) model. The test is asymptotic.
}

\usage{
AMUSEasymp(X, \ldots)

\method{AMUSEasymp}{default}(X, k, tau = 1, \ldots)
\method{AMUSEasymp}{ts}(X, \ldots)
\method{AMUSEasymp}{xts}(X, \ldots)
\method{AMUSEasymp}{zoo}(X, \ldots)

}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{k}{The number of latent series that are not white noise. Can be between \eqn{0} and \eqn{p-1}.}
  \item{tau}{The lag for the AMUSE autocovariance matrix.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
AMUSE standardizes \code{X} with \eqn{n} samples and computes the eigenedcomposition of the autocovariance matrix of the standardized data for a chosen lag \code{tau}, yielding a transformation \eqn{\bf W}{W} giving the latent variables as \eqn{{\bf S} = {\bf X} {\bf W}}{S = XW}. Assume, without loss of generality, that the latent components are ordered in decreasing order with respect to the squares of the corresponding eigenvalues of the autocovariance matrix. Under the null hypothesis the final \eqn{p - k} eigenvalues equal zero, \eqn{\lambda_{p-k} = \cdots = \lambda_{p} = 0}{\lambda_(p-k) = \ldots = \lambda_p = 0}, and their mean square \eqn{m} can be used as a test statistic in inference on the true number of latent white noise series. 

This function conducts the hypothesis test using the asymptotic null distribution of \eqn{m}, a chi-squared distribution with \eqn{(p - k)(p - k + 1)/2} degrees of freedom.
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
  \item{MU}{The location of the data which was subtracted before calculating AMUSE.}
  \item{tau}{The used lag.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\references{
Virta, J. and Nordhausen, K. (2019), \emph{Determining the Signal Dimension in Second Order Source Separation}. To appear in Statistica Sinica, <doi:10.5705/ss.202018.0347>.
}

\author{
Klaus Nordhausen, Joni Virta
}

\seealso{
\code{\link[JADE]{AMUSE}}, \code{\link[JADE]{SOBI}}, \code{\link{SOBIasymp}}
}

\examples{
  n <- 1000

  A <- matrix(rnorm(16), 4, 4)
  s1 <- arima.sim(list(ar = c(0.3, 0.6)), n)
  s2 <- arima.sim(list(ma = c(-0.3, 0.3)), n)
  s3 <- rnorm(n)
  s4 <- rnorm(n)

  S <- cbind(s1, s2, s3, s4)
  X <- S \%*\% t(A)

  asymp_res_1 <- AMUSEasymp(X, k = 1)
  asymp_res_1

  asymp_res_2 <- AMUSEasymp(X, k = 2)
  asymp_res_2

  # Plots of the estimated sources, the last two are white noise
  plot(asymp_res_2)
}

\keyword{ htest }
\keyword{ multivariate }
\keyword{ ts }