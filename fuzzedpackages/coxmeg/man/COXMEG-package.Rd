\name{coxmeg-package}
\alias{coxmeg-package}
\docType{package}
\title{
\packageTitle{coxmeg}
}
\description{
\packageDescription{coxmeg}
}
\details{

coxmeg is an R package for efficiently conducting GWAS of age-at-onset traits using a Cox mixed-effects model. coxmeg introduces a fast estimation algorithm for general sparse relatedness matrices including but not limited to block-diagonal pedigree-based matrices. coxmeg also introduces a fast and powerful score test for fully dense relatedness matrices, accounting for both population stratification and family structure. In addition, coxmeg can handle positive semidefinite relatedness matrices. Compared to coxme, coxmeg substantially improves the computational efficiency for estimating or testing genetic effects by using a variance component estimated once from a null model, and introducing fast algorithms, including inexact newton methods, preconditioned conjugate gradient methods and stochastic Lanczos quadrature.

}
\author{
\packageAuthor{coxmeg}

Maintainer: \packageMaintainer{coxmeg}
}
\references{
He. L, Kulminski, A.M., Genome-wide association analysis of age-at-onset traits using Cox mixed-effects models. 2019 
}

\keyword{ package }

\examples{
library(Matrix)
library(MASS)
library(coxmeg)

## simulate a block-diagonal relatedness matrix
tau_var <- 0.2
n_f <- 100
mat_list <- list()
size <- rep(10,n_f)
offd <- 0.5
for(i in 1:n_f)
{
   mat_list[[i]] <- matrix(offd,size[i],size[i])
   diag(mat_list[[i]]) <- 1
}
sigma <- as.matrix(bdiag(mat_list))
n <- nrow(sigma)

## simulate random effexts and outcomes
x <- mvrnorm(1, rep(0,n), tau_var*sigma)
myrates <- exp(x-1)
y <- rexp(n, rate = myrates)
cen <- rexp(n, rate = 0.02 )
ycen <- pmin(y, cen)
outcome <- cbind(ycen,as.numeric(y <= cen))

## fit the Cox mixed-effects model
re = coxmeg(outcome,sigma,type='bd',order=1,detap='diagonal')
re

}
