\name{tlrmvnmvt-package}
\alias{tlrmvnmvt}
\docType{package}
\title{
\packageTitle{tlrmvnmvt}
}
\description{
\packageDescription{tlrmvnmvt}
}
\details{
  Implementation of the classic Genz algorithm and a novel tile-low-rank algorithm for computing relatively high-dimensional multivariate normal and Student-t probabilities. For the Genz's algorithm (\code{GenzBretz}), we apply a univariate reordering preconditioner and for the tile-low-rank algorithms (\code{TLRQMC}), we apply a recursive block reordering preconditioner. The \code{GenzBretz} methods are different from their counterparts in the `mvtnorm` package in that the `tlrmvnmvt` package can accept any problem dimension and return the result in the log2 fashion, which is useful when the true probability is smaller than the machine precision. The \code{TLRQMC} algorithms can compute the probabilities up to tens of thousands of dimensions with the low-rank representation. However, this category of algorithms requires the existence of the low-rank property in the off-diagonal blocks of size \code{m}. The \code{zorder} function implements Morton's order in the 2D plane, which enhances the low-rank property of the produced covariance matrices.
  
  Package functions: pmvn, pmvt, zorder
}
\author{
\packageAuthor{tlrmvnmvt}

Maintainer: \packageMaintainer{tlrmvnmvt}
}
\references{
Cao, J., Genton, M. G., Keyes, D. E., & Turkiyyah, G. M. (2019). "Exploiting Low
Rank Covariance Structures for Computing High-Dimensional Normal and Student-
t Probabilities"
}
