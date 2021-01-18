\name{FBF_GS}


\alias{FBF_GS}


\title{Moment Fractional Bayes Factor Stochastic Search with Global Prior for Gaussian DAG Models}
\description{Estimate the edge inclusion probabilities for a Gaussian DAG with q nodes from observational data, using the moment fractional Bayes factor approach with global prior.
}


\usage{
FBF_GS(Corr, nobs, G_base, h, C, n_tot_mod, n_hpp)

}


\arguments{
  \item{Corr}{qxq correlation matrix.}
  \item{nobs}{Number of observations.}
  \item{G_base}{Base DAG.}
  \item{h}{Parameter prior.}
  \item{C}{Costant who keeps the probability of all local moves bounded away from 0 and 1.}
  \item{n_tot_mod}{Maximum number of different models which will be visited by the algorithm, for each equation.}
  \item{n_hpp}{Number of the highest posterior probability models which will be returned by the procedure.}
}


\value{An object of \code{class} \code{list} with:
 \describe{
    \item{\code{M_q}}{
      Matrix (qxq) with the estimated edge inclusion probabilities.
    }
    \item{\code{M_G}}{
      Matrix (n*n_hpp)xq with the n_hpp highest posterior probability models returned by the procedure.
    }
    \item{\code{M_P}}{
      Vector (n_hpp) with the n_hpp posterior probabilities of the models in M_G.
    }
  }
}


\references{
  D. Altomare, G. Consonni and L. La Rocca (2012). Objective Bayesian search of Gaussian directed acyclic graphical models for ordered variables with non-local priors. \emph{Article submitted to Biometric Methodology}.
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

Res_search=FBF_GS(Corr, nobs, matrix(0,q,q), 1, 0.01, 1000, 10)
M_q=Res_search$M_q
M_G=Res_search$M_G
M_P=Res_search$M_P

G_med=M_q
G_med[M_q>=0.5]=1
G_med[M_q<0.5]=0 #median probability DAG

G_high=M_G[1:q,1:q] #Highest Posterior Probability DAG (HPP)
pp_high=M_P[1] #Posterior Probability of the HPP

#Structural Hamming Distance between the true DAG and the median probability DAG
sum(sum(abs(G_med-Gt)))
#Structural Hamming Distance between the true DAG and the highest probability DAG 
sum(sum(abs(G_high-Gt)))

## End(Not run) 


}


\keyword{multivariate}
\keyword{models}
\keyword{dag}
\keyword{stochastic search}
