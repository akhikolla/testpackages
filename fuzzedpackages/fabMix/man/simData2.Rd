\name{simData2}
\alias{simData2}
\title{
Synthetic data generator 2
}
\description{
Simulate data from a multivariate normal mixture using a mixture of factor analyzers mechanism. 
}
\usage{
simData2(sameSigma,  p, q, K.true, n, loading_means, loading_sd, sINV_values)
}
\arguments{
  \item{sameSigma}{
	Logical. 
}
  \item{p}{
	The dimension of the multivariate normal distribution (\eqn{p > 1}). 
}
  \item{q}{
	Number of factors. It should be strictly smaller than p.
}
  \item{K.true}{
	The number of mixture components (clusters).
}
  \item{n}{
	Sample size.
}
  \item{loading_means}{
	A vector which contains the means of blocks of factor loadings.

	Default: \code{loading_means = c(-30,-20,-10,10, 20, 30)}.
}
  \item{loading_sd}{
	A vector which contains the standard deviations of blocks of factor loadings. 

	Default: \code{loading_sd <- rep(2, length(loading_means))}.
}
  \item{sINV_values}{
	A vector which contains the values of the diagonal of the (common) inverse covariance matrix, if \code{sigmaTrue = TRUE}. An \eqn{K\times p} matrix which contains the values of the diagonal of the inverse covariance matrix per component, if \code{sigmaTrue = FALSE}. 

	Default: \code{ sINV_values = rgamma(p, shape = 1, rate = 1)}.
}
}
\value{
	A list with the following entries:
	\item{data}{
                \eqn{n\times p} array containing the simulated data.
        }
	\item{class}{
                \eqn{n}-dimensional vector containing the class of each observation.
        }
	\item{factorLoadings}{
                \eqn{K.true\times p \times q}-array containing the factor loadings \eqn{\Lambda_{krj}} per cluster \eqn{k}, feature \eqn{r} and factor \eqn{j}, where \eqn{k=1,\ldots,K}; \eqn{r=1,\ldots,p}; \eqn{j=1,\ldots,q}.
        }
	\item{means}{
		\eqn{K.true\times p} matrix containing the marginal means \eqn{\mu_{kr}}, \eqn{k=1,\ldots,K}; \eqn{r=1,\ldots,p}.
	}
	\item{variance}{
		\eqn{p\times p} diagonal matrix containing the variance of errors \eqn{\sigma_{rr}}, \eqn{r=1,\ldots,p}. Note that the same variance of errors is assumed for each cluster.
	}
	\item{factors}{
		\eqn{n\times q} matrix containing the simulated factor values.
	}
	\item{weights}{
		\eqn{K.true}-dimensional vector containing the weight of each cluster.
	}

}
\author{
Panagiotis Papastamoulis
}
\note{
The marginal variance for cluster \eqn{k} is equal to \eqn{\Lambda_k\Lambda_k^{T} + \Sigma}.}
\examples{
library('fabMix')

n = 8                # sample size
p = 5                # number of variables
q = 2                # number of factors
K = 2                # true number of clusters

sINV_diag = 1/((1:p))    # diagonal of inverse variance of errors
set.seed(100)
syntheticDataset <- simData2(K.true = K, n = n, q = q, p = p, 
                        sINV_values = sINV_diag)
summary(syntheticDataset)
}

