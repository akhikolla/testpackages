\name{gmd}
\title{Gini Mean Difference}
\alias{gmd}
\description{
  Computes Gini mean difference of x, where alpha is an exponent on the Euclidean distance and return the Gini mean difference. The default value for alpha is 1. 
}
\usage{
  gmd(x, alpha)
}
\arguments{
  \item{x}{data}
  \item{alpha}{exponent on Euclidean distance, in (0,2)}
}
\details{
  \code{gmd} compute Gini mean difference of data.
  It is a self-contained R function dealing with both univariate and multivariate data. 
  
  The samples must not contain missing values.  \code{alpha} if missing by default is 1, otherwise it is exponent on the Euclidean distance.
  
  Gini mean difference (GMD) was originally introduced as an alternative measure of variability to the usual standard deviation (\cite{Gini14}, \cite{Yitzhaki13}). Let \eqn{X} and \eqn{X^\prime} be independent random variables from a univariate distribution \eqn{F} with finite first moment in \eqn{R}.  The GMD of \eqn{F} is
    
    \eqn{\Delta=\Delta(X)=\Delta(F)=E|X-X^{\prime}|,}
  
  the expected distance between two independent random variables. If the sample data \eqn{\mathbf x=\{x_1,x_2,...,x_n\}} is available, the sample Gini mean difference is calculated by 
  
  \eqn{\hat{\Delta} = {n \choose 2}^{-1} \sum_{1\leq i<j\leq n} | x_i - x_j| = {n \choose 2}^{-1} \sum_{i=1}^n (2i-n-1) x_{(i)},}

where \eqn{x_{(1)} \leq x_{(2)} \leq \cdots \leq x_{(n)}} are the order statistics of \eqn{\mathbf x} (\cite{Schechtman87}). The computation complexity for univariate Gini Mean difference is \eqn{O(n \log n)}.

Gini mean difference has been generalized for multivariate distributions (\cite{Koshvoy97})  That is, the Gini mean difference of a distribution F in \eqn{\mathbf{R}^d} is  \eqn{\Delta =E \|\mathbf X -\mathbf X ^\prime\|,}
or even more generally for some \eqn{\alpha \in (0,2)},

\eqn{\Delta(\alpha) = E \|\mathbf X-\mathbf X^\prime\|^{\alpha}}, 

where \eqn{\| \mathbf x \|} is the Euclidean norm. The sample Gini mean difference is computed by 
  
  \eqn{\hat{\Delta(\alpha)} = {n \choose 2}^{-1} \sum_{1\leq i<j\leq n} \| x_i - x_j\|^{\alpha}.}

Its computation complexity is \eqn{O(n^2)}.


}
\value{
  \code{gmd} returns the sample Gini mean distance. 
}
\seealso{
  \code{\link{RcppGmd}}  \code{\link{gCov}}  \code{\link{gCor}}  
}

\references{
  
  Gini, C. (1914). Sulla misura della concentrazione e della variabilita dei caratteri. Atti del Reale Istituto Veneto di Scienze, Lettere ed Aeti, 62, 1203-1248. English Translation: On the measurement of concentration and variability of characters (2005). Metron, LXIII(1), 3-38.
  
  
  Koshevoy, G. and Mosler, K. (1997). Multivariate Gini indices. Journal of Multivariate Analysis, 60, 252-276.
  
  Schechtman, E. and Yitzhaki, S. (1987). A measure of association based on Gini's mean difference. Communication in Statistics-Theory and Methods, 16 (1), 207-231.
  
  Yitzhaki, S. and Schechtman, E. (2013). The Gini Methodology, Springer, New York. 
} 

\examples{
  n = 100
  x <- runif(n)
  
  t0 = proc.time()
  gmd(x, alpha=1)
  proc.time()- t0
  
  t1 = proc.time()
  gmd(x, alpha=0.5)    
  proc.time()- t1
  
  x <- matrix(runif(n), n/2, 2)
  gmd(x,alpha=1)
  
}

\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
