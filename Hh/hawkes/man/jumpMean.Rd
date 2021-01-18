\name{jumpMean}
\alias{jumpMean}
\title{
Mean of Hawkes process jumps.
}
\description{
The function returns the theoretical mean of the number of jumps of a Hawkes process on a time interval of length tau.}
\usage{
jumpMean(lambda0, alpha, beta, tau)
}

\arguments{
  \item{lambda0}{Vector of initial intensity, a scalar in the monovariate case.}
  \item{alpha}{Matrix of excitation, a scalar in the monovariate case. Excitation values are all positive.}
  \item{beta}{Vector of betas, a scalar in the monovariate case.}
  \item{tau}{Time interval length.}
}
\details{
 Notice that in the scalar case, one must have beta>alpha for the process to be stable, and in the multivariate case, the matrix (diag(beta)-alpha) must have eigen values with strictly positive real parts for the process to be stable.
}
\value{
Returns a vector containing the mean number of jumps of every process component.
}
\references{
Jose Da Fonseca and Riadh Zaatour
Hawkes Process : Fast Calibration, Application to Trade Clustering and Diffusive Limit.
\emph{Journal of Futures Markets}, Volume 34, Issue 6, pages 497-606, June 2014.

Jose Da Fonseca and Riadh Zaatour
Clustering and Mean Reversion in Hawkes Microstructure Models.

}

\examples{
#One dimensional Hawkes process
lambda0<-0.02
alpha<-0.05
beta<-0.06
tau<-60#one minute
h<-jumpMean(lambda0,alpha,beta,tau)

#Multivariate Hawkes process
lambda0<-c(0.02,0.02)
alpha<-matrix(c(0.05,0,0,0.05),byrow=TRUE,nrow=2)
beta<-c(0.06,0.06)
tau<-60#one minute
h<-jumpMean(lambda0,alpha,beta,tau)
}
