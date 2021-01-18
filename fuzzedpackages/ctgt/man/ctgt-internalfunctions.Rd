\name{internalfunctions}
\alias{actgt_it}
\alias{tacmax}
\alias{tactrue}

\title{'Internal Functions (ctgt)'}
\description{Internal functions of ctgt.}
\usage{
## iterative shortcut with branch and bound
actgt_it(y,Tmatrix, Cmatrix,fxs, sxs,Tf,Lamf,Cf,Ts,Lams,Cs,count=1,maxIt=1,a = 0.05)
  
## to check whether tmin is above cmax
tacmax(tmins,levels,tw, cf,lf,ls,alp )
  
## to check whether tmin is above ctrue
tactrue(tmins,hyxs,cfull,Wmatrix,alp )


}

\arguments{
  \item{y}{The response vector (numeric vector).}
  
  \item{Tmatrix}{The matrix used to calculate the test statistics.}
  
  \item{Cmatrix}{The matrix used to calculate the critical values.}
  
  \item{fxs}{The name vector of upper model (character vector). }
  
  \item{sxs}{The name vector of lower model (character vector). }
  
  \item{Tf,Lamf,Cf}{Test statistic, eigenvalues, and critical value of fxs.}
  
  \item{Ts,Lams,Cs}{Test statistic, eigenvalues, and critical value of sxs.}
  
  \item{count}{count the branches, default is 1.}
  
  \item{maxIt}{maximal number of branches chosen by user, default is 1.}
  
  \item{a, alp}{alpha level.}
  
  \item{tmins}{Minimum test statistics. }
  
  \item{levels}{levels}
  
  \item{tw, cf,lf,ls}{sorted weights, critical values and level for fxs and sxs. }
  
  \item{hyxs}{The name vector of the covariates of interest (character vector).}
  
  \item{cfull}{ critical value of full model.}
  
  \item{Wmatrix}{matrix to calculate majorizing vector.}
  
}


\author{
Ningning Xu

Maintainer: Ningning Xu <n.xu@lumc.nl; xu15263142750@gmail.com>
}

