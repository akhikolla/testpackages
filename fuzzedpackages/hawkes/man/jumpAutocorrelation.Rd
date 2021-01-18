\name{jumpAutocorrelation}
\alias{jumpAutocorrelation}
\title{
Autocorrelation of Hawkes process jumps on nonoverlapping time intervals with lag.
}
\description{
The function returns the theoretical autocorrelation of the number of jumps of a Hawkes process on nonoverlapping time intervals with lag}
\usage{
jumpAutocorrelation(lambda0, alpha, beta, tau,lag)
}

\arguments{
  \item{lambda0}{Vector of initial intensity, a scalar in the monovariate case.}
  \item{alpha}{Matrix of excitation, a scalar in the monovariate case. Excitation values are all positive.}
  \item{beta}{Vector of betas, a scalar in the monovariate case.}
  \item{tau}{Time interval length.}
  \item{lag}{Time lag.}
}
\details{
 Notice that in the scalar case, one must have beta>alpha for the process to be stable, and in the multivariate case, the matrix (diag(beta)-alpha) must have eigen values with strictly positive real parts for the process to be stable.
}
\value{
Returns a matrix containing the autocorrelation of the number of jumps of process components.
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
lag<-0#adjacent non overlappingintervals
h<-jumpAutocorrelation(lambda0,alpha,beta,tau,lag)

#Multivariate Hawkes process
lambda0<-c(0.02,0.02)
alpha<-matrix(c(0.05,0,0,0.05),byrow=TRUE,nrow=2)
beta<-c(0.06,0.06)
tau<-60#one minute
lag<-0#adjacent non overlappingintervals
h<-jumpAutocorrelation(lambda0,alpha,beta,tau,lag)
}
