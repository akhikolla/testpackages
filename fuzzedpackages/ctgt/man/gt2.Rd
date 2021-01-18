\name{non-standardized globaltest}
\alias{gt2}

\title{'Non-standardized globaltest'}

\description{This is the sencond version of the globaltest, the non-standardized globaltest}

\usage{
  ## a powerful variant of globaltest
  gt2 (y, X, hyps, alpha = 0.05)

}

\arguments{
  \item{y}{The response vector (numeric vector).}
  \item{X}{The full design matrix, whose columns are named by the covariates.}
  \item{hyps}{The name vector of the covariates in the pathway of interest (character vector).}
  \item{alpha}{The type I error rate allowed. The default is 0.05.}
}

\value{Returns the p-value, the observed and expected test statistics and the number of covariates. }

\author{
Ningning Xu

Maintainer: Ningning Xu <n.xu@lumc.nl; xu15263142750@gmail.com>
}
\references{
   Ningning Xu, Aldo solari, Jelle Goeman, Clsoed testing with global test, with applications on metabolomics data, arXiv:2001.01541, https://arxiv.org/abs/2001.01541
}

\examples{
  #Generate the design matrix and response vector for logistic regression models
  n= 100
  m = 5
  X = matrix(0, n, m,byrow = TRUE )
  for ( i in 1:n){
    set.seed(1234+i)
    X[i,] =  as.vector(arima.sim(model = list(order = c(1, 0, 0), ar = 0.2), n = m) )
  }

  y = rbinom(n,1,0.6)
  X[which(y==1),1:3] = X[which(y==1),1:3] + 0.8


  xs = paste("x",seq(1,m,1),sep="") 
  colnames(X) = xs

  hyps=xs[1]

  #The raw p-values of globaltest
  gt2(y,X,hyps)
  #p-value Statistic  Expected      #Cov 
  #7.64e-03  2.30e+02  1.24e+02  1.00e+00 
}