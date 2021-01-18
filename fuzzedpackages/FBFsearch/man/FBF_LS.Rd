\name{FBF_LS}


\alias{FBF_LS}


\title{Moment Fractional Bayes Factor Stochastic Search with Local Prior for DAG Models}
\description{Estimate the edge inclusion probabilities for a directed acyclic graph (DAG) from observational data, using the moment fractional Bayes factor approach with local prior.
}


\usage{
FBF_LS(Corr, nobs, G_base, h, C, n_tot_mod)

}


\arguments{
  \item{Corr}{qxq correlation matrix.}
  \item{nobs}{Number of observations.}
  \item{G_base}{Base DAG.}
  \item{h}{Parameter prior.}
  \item{C}{Costant who keeps the probability of all local moves bounded away from 0 and 1.}
  \item{n_tot_mod}{Maximum number of different models which will be visited by the algorithm, for each equation.}
}


\value{An object of \code{class} \code{matrix} with the estimated edge inclusion probabilities.}


\references{
  D. Altomare, G. Consonni and L. LaRocca (2012).Objective Bayesian search of Gaussian directed acyclic graphical models for ordered variables with non-local priors.\emph{Article submitted to Biometric Methodology}.
}


\author{
  Davide Altomare (\email{davide.altomare@gmail.com}).
}


\examples{

## Not run:

data(SimDag6) 

Corr=dataSim6$SimCorr[[1]]
nobs=50
q=ncol(Corr)
Gt=dataSim6$TDag

M_q=FBF_LS(Corr, nobs, matrix(0,q,q), 0, 0.01, 1000)

G_med=M_q
G_med[M_q>=0.5]=1
G_med[M_q<0.5]=0 #median probability DAG

#Structural Hamming Distance between the true DAG and the median probability DAG
sum(sum(abs(G_med-Gt))) 

## End(Not run)

}


\keyword{multivariate}
\keyword{models}
\keyword{dag}
\keyword{stochastic search}
