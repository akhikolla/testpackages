\name{SOBIladle}
\alias{SOBIladle}
\alias{SOBIladle.default}
\alias{SOBIladle.ts}
\alias{SOBIladle.xts}
\alias{SOBIladle.zoo}
\title{
Ladle Estimator to Estimate the Number of White Noise Components in SOS with SOBI
}

\description{
The ladle estimator uses the joint diagonalization "eigenvalues" and "eigenvectors" of several autocovariance matrices to estimate the number of white noise components in second-order source separation (SOS).
}

\usage{
SOBIladle(X, \ldots)

\method{SOBIladle}{default}(X, tau = 1:12, l = 20, sim = c("geom", "fixed"), n.boot = 200,
ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1),
maxiter = 1000, eps = 1e-06, \ldots)
\method{SOBIladle}{ts}(X, \ldots)
\method{SOBIladle}{xts}(X, \ldots)
\method{SOBIladle}{zoo}(X, \ldots)

}

\arguments{
  \item{X}{A numeric matrix or a multivariate time series object of class \code{\link{ts}}, \code{\link[xts]{xts}} or \code{\link[zoo]{zoo}}. Missing values are not allowed.}
  \item{tau}{The lags for the SOBI autocovariance matrices.}
  \item{l}{If \code{sim = "geom"} then \code{l} is the success probability of the geometric distribution from where the bootstrap block lengths for the stationary bootstrap are drawn. If \code{sim = "fixed"} then \code{l} is the fixed block length for the fixed block bootstrap.}
  \item{sim}{If \code{"geom"} (default) then the stationary bootstrap is used. If \code{"fixed"} then the fixed block bootstrap is used.}
  \item{n.boot}{The number of bootstrapping samples. See \code{\link[boot]{tsboot}} for details.}
  \item{ncomp}{The number of components among which the ladle estimator is to be searched. Must be between \code{0} and \code{ncol(X)-1}. The default follows the recommendation of Luo and Li (2016).}
  \item{maxiter}{Maximum number of iterations.}
  \item{eps}{Convergence tolerance.}
  \item{\ldots}{Further arguments to be passed to or from methods.}
}

\details{
SOBI standardizes \code{X} with \eqn{n} samples and jointly diagonalizes the autocovariance matrices of the standardized data for a chosen set of lags \code{tau}, yielding a transformation \eqn{\bf W}{W} giving the latent variables as \eqn{{\bf S} = {\bf X} {\bf W}}{S = XW}. Assume, without loss of generality, that the latent components are ordered in decreasing order with respect to the sums of squares of the corresponding "eigenvalues" produced by the joint diagonalization.
Under the assumption that we have \eqn{k} non-white-noise components, the final \eqn{p - k} "eigenvalues" of each of the autocovariance matrices equal zero, \eqn{\lambda^\tau_{p-k} = \cdots = \lambda^\tau_{p} = 0}{\lambda^tau_(p-k) = \ldots = \lambda^tau_p = 0}.

The change point from non-zero eigenvalues to zero eigenvalues is visible in the joint diagonalization "eigenvectors" of the autocovariance matrices as an increase in their boostrap variablity. Similarly, before the change point, the squared eigenvalues decrease in magnitude and afterwards they stay constant. The ladle estimate combines the scaled eigenvector bootstrap variability with the scaled eigenvalues to estimate the number of non-white-noise components. The estimate is the value of \eqn{k = 0, \ldots ,} \code{ncomp} where the combined measure achieves its minimum value. 

This function uses for the joint diagonalization the function \code{\link[JADE:rjd]{frjd.int}}, which does not fail in case of failed convergence but returns the estimate from the final step.
}

\value{
A list of class \code{ladle} containing:
  \item{method}{The string SOBI.}
  \item{k}{The estimated number of non-white-noise components.}
  \item{fn}{The vector giving the measures of variation of the eigenvectors using the bootstrapped eigenvectors for the different number of components.}
  \item{phin}{Normalized sums of squared eigenvalues of the SOBI matrices.}
  \item{data.name}{The name of the data for which the ladle estimate was computed.}
  \item{gn}{The main criterion for the ladle estimate - the sum of \code{fn} and \code{phin}. \code{k} is the value where \code{gn} takes its minimum.}
  \item{lambda}{The sums of squared eigenvalues of the SOBI matrices.}
  \item{W}{The transformation matrix to the source components. Also known as the unmixing matrix.}
  \item{S}{Multivariate time series with the centered source components.}
  \item{MU}{The location of the data which was subtracted before calculating the source components.}
  \item{sim}{The used boostrapping technique, either \code{"geom"} or \code{"fixed"}.}
  \item{tau}{The used set of lags for the SOBI autocovariance matrices.}
}

\references{
Nordhausen, K. and Virta, J.(2018), \emph{Ladle Estimator for Time Series Signal Dimension}. In 2018 IEEE Statistical Signal Processing Workshop (SSP), pp. 428--432, <doi:10.1109/SSP.2018.8450695>.
  
Luo, W. and Li, B. (2016), \emph{Combining Eigenvalues and Variation of Eigenvectors for Order Determination}, Biometrika, 103. 875--887. <doi:10.1093/biomet/asw051>
}

\author{
Klaus Nordhausen, Joni Virta
}
\seealso{
\code{\link[JADE]{AMUSE}}, \code{\link[JADE]{SOBI}}, \code{\link{AMUSEladle}}, \code{\link[JADE:rjd]{frjd.int}}
}

\examples{
  n <- 1000
  
  s1 <- arima.sim(n = n, list(ar = 0.6, ma = c(0, -0.4)))
  s2 <- arima.sim(n = n, list(ar = c(0, 0.1,0.3), ma = c(0.2, 0.4)))
  s3 <- arima.sim(n = n, list(ar = c(0, 0.8)))
  Snoise <- matrix(rnorm(5*n), ncol = 5)
  S <- cbind(s1, s2, s3, Snoise)

  A <- matrix(rnorm(64), 8, 8)
  X <- S \%*\% t(A)
  
  ladle_SOBI <- SOBIladle(X, l = 20, sim = "geom")

  # The estimated number of non-white-noise components
  summary(ladle_SOBI)
  
  # The ladle plot
  ladleplot(ladle_SOBI)
  
  # Time series plots of the estimated components
  plot(ladle_SOBI)
  
  # Note that AMUSEladle with lag 1 does not work due to the lack of short range dependencies
  ladle_AMUSE <- AMUSEladle(X)

  summary(ladle_AMUSE)
  ladleplot(ladle_AMUSE)
  
  # xts series as input
  library(xts)
  data(sample_matrix)
  X2 <- as.xts(sample_matrix)
  ladle_SOBI_xts <- SOBIladle(X2, l = 20, sim = "geom")
  plot(ladle_SOBI_xts, multi.panel = TRUE)

  # zoo series as input
  X3 <- as.zoo(X)
  ladle_SOBI_zoo <- SOBIladle(X3, l = 20, sim = "geom")
  plot(ladle_SOBI_zoo)
}

\keyword{ multivariate }
\keyword{ ts }