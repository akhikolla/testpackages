\name{symm}
\alias{symm}
\title{
Symmetry Statistics for all Interactions
}
\description{
\code{symm} returns all the symmetry statistics up to depth \code{d} in marginal binary expansions for the tests \code{BET} and \code{BETs}.
}
\usage{
symm(X, d, unif.margin = FALSE)
}
\arguments{
  \item{X}{
a matrix to be tested.
}
  \item{d}{
depth of the marginal binary expansions.
}
  \item{unif.margin}{
  logicals. If \code{TRUE} the data has been uniformed based on empirical cumulative distribution function. Default to be \code{FALSE} and the function uniforms the data.
  }
}

\value{
\item{SymmetryStatistics}{a dataframe with \eqn{(p+2)} columns, where \eqn{p} is the number of columns of \code{X}. The first column gives the binary index for all variables, the next \eqn{p} columns displays all the interactions of respective variables, the last column of \code{SymmetryStatistics} gives the respective symmetry statistic.}
}

\examples{
v <- runif(128, -pi, pi)
X1 <- cos(v) + 2.5 * rnorm(128, 0, 1/20)
X2 <- sin(v) + 2.5 * rnorm(128, 0, 1/20)
symm(cbind(X1, X2), 3)
}


