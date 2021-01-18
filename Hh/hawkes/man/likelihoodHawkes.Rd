\name{likelihoodHawkes}
\alias{likelihoodHawkes}

\title{
Compute the likelihood function of a hawkes process
}
\description{
Compute the likelihood function of a hawkes process for the given parameter and given the jump times vector (or list of vectors in the multivariate case), and until a time horizon.
}
\usage{
likelihoodHawkes(lambda0, alpha, beta, history)
}

\arguments{
  \item{lambda0}{Vector of initial intensity, a scalar in the monovariate case.}
  \item{alpha}{Matrix of excitation, a scalar in the monovariate case. Excitation values are all positive.}
  \item{beta}{Vector of betas, a scalar in the monovariate case.}
  \item{history}{Jump times vector (or list of vectors in the multivariate case).}
}

\value{
Returns the opposite of the likelihood.

}
\references{
Y. Ogata. (1981)
On Lewis simulation method for point processes.
\emph{IEEE Transactions on Information Theory}, \bold{31}
}

\examples{
#One dimensional Hawkes process
lambda0<-0.2
alpha<-0.5
beta<-0.7
history<-simulateHawkes(lambda0,alpha,beta,3600)
l<-likelihoodHawkes(lambda0,alpha,beta,history[[1]])

#Multivariate Hawkes process
lambda0<-c(0.2,0.2)
alpha<-matrix(c(0.5,0,0,0.5),byrow=TRUE,nrow=2)
beta<-c(0.7,0.7)
history<-simulateHawkes(lambda0,alpha,beta,3600)
l<-likelihoodHawkes(lambda0,alpha,beta,history)
}
