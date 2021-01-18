\name{sl0fused_c}
\alias{sl0fused_c}
\title{Sparse L0 Fused Regression}
\description{
  This is a cpp function used for R function sl0fused.
}

\usage{
 sl0fused_c(y, T0, T02, max_steps)
}
\arguments{
  \item{y}{Response sequence to be fitted.}
  \item{T0}{Number of knots in the piecewise constant(breaks in the derivative).}
  \item{T02}{Number of non-zero values in fitted coefficient.}
  \item{max_steps}{The maximum number of steps for the algorithm to take before termination. }
}
\value{
  \item{beta}{Fitted value.}
 \item{u}{Dual coefficient.}
 \item{z}{Primal coefficient.}
 \item{d}{Dual coefficient.}
}


\author{Canhong Wen, Xueqin Wang, Yanhe Shen, Aijun Zhang}
\references{
 Wen,C., Wang, X., Shen, Y., and Zhang, A. (2017). "L0 trend filtering", technical report.
}
