\name{l0tf_c}
\alias{l0tf_c}
\title{L0 Trend Filtering}
\description{
This is a cpp function used for R function l0tf.
}

\usage{
 l0tf_c(y, k0, T0, max_steps)
}
\arguments{
  \item{y}{Numeric vector of inputs.}
  \item{k0}{An integer specifying the desired order of the piecewise polyomial produced by the solution of the trend filtering problem. }
  \item{T0}{Number of knots in the piecewise polynomial(breaks in the (k+1)-st derivative).}
  \item{max_steps}{The maximum number of steps for the algorithm to take before termination. }
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
