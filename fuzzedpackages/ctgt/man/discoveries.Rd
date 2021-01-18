\name{true discoveries}
\alias{discoveries}

\title{True discoveries}

\description{To count the number of true discoveries within a given pathway or feature set of interest.}

\usage{
  discoveries (y, X, xs, hyps, maxit = 0, alpha = 0.05)
}

\arguments{
  \item{y}{The response vector (numeric vector).}
  \item{X}{The full design matrix, whose columns are named by the covariates.}
  \item{xs}{The name vector of all covariates (character vector). }
  \item{hyps}{The name vector of the covariates in the pathway of interest (character vector).}
  \item{maxit}{An optional integer to denote the maximal interations for branch and bound method. The default value 0 means the single-step shortcut without branch and bound method. Note that larger value is more time-consuming. }
  \item{alpha}{The type I error rate allowed. The default is 0.05.}
  
  
}

\value{Returns a non-negative interger. }

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

  # For standarized data
  X = scale(X,center = FALSE,scale = TRUE)/sqrt(n-1)
  interest = xs

  discoveries(y,X,xs,hyps = interest)
  #2
  discoveries(y,X,xs,hyps = interest,maxit=10)
  #2
}