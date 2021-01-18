
\name{mvnintGHK}
\alias{mvnintGHK}
\title{Computing Multivariate Normal Rectangle Probability}
\usage{
mvnintGHK(mean, sigma, lower, upper, nrep = 5000, log = TRUE,
  reorder = TRUE)
}
\arguments{
\item{mean}{the numeric vector of mean of length \code{n}.}

\item{sigma}{the covariance matrix of dimension \code{n}.}

\item{lower}{the numeric vector of lower limits of length \code{n}.}

\item{upper}{the numeric vector of upper limits of length \code{n}.}

\item{nrep}{a positive integer of Monte Carlo size.}

\item{log}{if TRUE then return the log of the probability. If FALSE return the probability.}

\item{reorder}{if TRUE then variable reordering algorithm is applied. If FALSE then original ordering is used.}
}
\value{
A list of the following two components:

\item{value}{the value of the integral. If \code{log = TRUE} then output the log of the integral. }

\item{error}{ the Monte Carlo standard deviation.}

}
\description{
Computes the multivariate normal rectangle probability for arbitrary
             limits and covariance matrices using (reordered) sequential importance sampling.
}
\details{
This program implemented the Geweke-Hajivassiliou-Keane simulator of computing the
         multivariate normal rectangle probability. For more details see Keane (1994).
         Also a variable reordering algorithm in Gibson, etal (1994) was implemented.

         Note that both \code{-Inf} and \code{Inf} may be specified in \code{lower} and \code{upper}.
}
\examples{
mvnintGHK(mean = rep(0, 51), sigma =  diag(0.2, 51) + matrix(0.8, 51, 51),
          lower = rep(-2,51), upper = rep(2,51), nrep = 10000)
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\references{
Gibson GJ., Glasbey CA. and Elston DA. (1994) Monte Carlo evaluation of multivariate normal integrals and
 sensitivity to variate ordering. \emph{Advances in Numerical Methods and Applications,
 World Scientific Publishing, River Edge}.

 Keane, M. (1994) A computationally practical simulation estimator for panel data.
 \emph{Econometrica, 62:95-116}.
}
\seealso{
\code{\link[mvtnorm]{pmvnorm}}
}

\keyword{Estimation}
