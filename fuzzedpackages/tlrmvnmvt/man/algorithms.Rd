\name{algorithm}
\alias{GenzBretz}
\alias{TLRQMC}
\title{Parameters for the Quasi-Monte Carlo sampling}
\description{
  The two functions return objects containing the parameters for the dense-matrix based Quasi-Monte Carlo method and the tile-low-rank-matrix based Quasi-Monte Carlo method, respectively.
}
\usage{
  GenzBretz(N = 499) 
  TLRQMC(N = 499, m = 64, epsl = 1e-4)
}
\arguments{
  \item{N}{an integer, specifying the number of per-batch Monte Carlo samples. The 
    total number of Monte Carlo samples is 20 X N}
  \item{m}{an integer, specifying the block size for the tile-low-rank methods}
  \item{epsl}{numeric value, specifying the truncation level for the tile-low-rank methods}
}
\value{
  Return an object of the class "GenzBretz" or "TLRQMC", which is used as the parameter of the \code{\link{pmvn}} and the \code{\link{pmvt}} functions.
}
\author{Jian Cao, Marc Genton, David Keyes, George Turkiyyah}

