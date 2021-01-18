\name{loglikfuncmmmkmv}

\alias{loglikfuncmmmkmv}


\title{Calculate the loglikeligood for a general mixed model}
\description{
Calculate the loglikeligood for the LMM model expressed as \deqn{Y=XB+\sum_{j=1}^k Z_j G_j+WE,} where \eqn{Y} is the \eqn{n \times d} response variable, \eqn{X} is the \eqn{n \times q} design matrix of \eqn{q \times d} the fixed effects \eqn{B}, \eqn{Z_j} for \eqn{j=1,2,\ldots,k}  (\eqn{k\geq 1}) are the \eqn{n \times q_j} design matrices of the \eqn{q_j \times d} random effects \eqn{G_j,}  and \eqn{W} is the \eqn{n \times t} design matrix of \eqn{t \times d} residual effecs \eqn{E}. The random effects and the residual are independently distributed, and have matrix variate distributions (\eqn{G_j\sim N_{q_j \times d}(0_{q_j \times d}, K_j,\Sigma_j)} for \eqn{j=1,2,\ldots,k} and \eqn{E\sim N_{t \times d}(0_{t \times d}, R,\Sigma_E)}). 
}
\usage{
loglikfuncmmmkmv(Y,X,Zlist, Klist, sigmahatlist, B,W,R )
}
\arguments{
  \item{Y}{a numeric vector for the parameters (a mapping of the original parameters)}
  \item{X}{a numeric matrix (see examples and details)}
  \item{Zlist}{a numeric matrix (see examples and details)}
  \item{Klist}{a numeric matrix (see examples and details)}
 \item{sigmahatlist}{a numeric matrix (see examples and details)}
 \item{B}{bla bla}
 \item{W}{bla bal}
 \item{R}{bla bla}
}
\value{
bla bla
}
\details{

bla bla
}

\examples{
\dontrun{
library(SAMM)
n=100
nsample=80
rhotrans=5
ar1cov_cppforR(c(rhotrans),matrix(5))
rho=(2/pi)*atan(rhotrans)
rho
tan((pi/2)*(rho))

M1<-matrix(rbinom(n*300, 2, .2)-1, nrow=n)
K1<-relmatcov_cppforR(c(.01), M1)

M2<-matrix(rbinom(n*300, 2, .2)-1, nrow=n)
K2<-relmatcov_cppforR(c(0.03), M2)
W=(diag(5)[sample(1:5,n, replace=TRUE),])
covY<-3*K1+5*K2+10*(W\%*\%ar1cov_cppforR(c(rhotrans),matrix(5))\%*\%t(W))
K1[1:5,1:5]
dim(W)
dim(ar1cov_cppforR(c(6),matrix(5)))
Y<-10+crossprod(chol(covY),rnorm(n))


#training set
Trainset<-sample(1:n,nsample)
ytrain=Y[Trainset]
Xtrain=matrix(rep(1, n)[Trainset], ncol=1)
Ztrain=diag(n)[Trainset,]
Wtrain=W[Trainset,]

samout<-SAMM(Y=matrix(ytrain,ncol=1),X=Xtrain,
Zlist=list(Ztrain, Ztrain), Klist=list(K1,K2),
lambda=0, W=Wtrain,R=list("ar1",c(0),matrix(5,1,1)),
Siglist=list("","",""), corfunc=c(F,F,T), corfuncfixed=c(F,F,F),
sigfunc=c(F,F,F),mmalg="dermm_reml2", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)
samout$corfuncparamslist[[3]]
rhohat=(2/pi)*atan(samout$corfuncparamslist[[3]])
rhohat
ar1cov_cppforR(c(samout$corfuncparamslist[[3]]),matrix(5,1,1))
}
}

\author{Deniz Akdemir}

