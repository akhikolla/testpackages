\name{l0gen_c}
\alias{l0gen_c}
\title{L0 Generalized Regression}
\description{
  This is a cpp function used for R function l0gen.
}

\usage{
 l0gen_c(y, D, T0, max_steps, ddinv)
}
\arguments{
  \item{y}{Response sequence to be filtered.}
  \item{D}{Penalty matrix on coeffient beta.}
  \item{T0}{Number of knots in the penalized coefficient(breaks in the \code{D*beta}), same as s.}
  \item{max_steps}{The maximum number of steps for the algorithm to take before termination. Same as K.max.}
  \item{ddinv}{The inverse matrix of \code{D*t(D)}, could be NULL input.}
}
\value{
  \item{beta}{Fitted value.}
 \item{u}{Dual coefficient.}
 \item{z}{Primal coefficient.}
}

\author{Canhong Wen, Xueqin Wang, Yanhe Shen, Aijun Zhang}
\references{
 Wen,C., Wang, X., Shen, Y., and Zhang, A. (2017). "L0 trend filtering", technical report.
}


