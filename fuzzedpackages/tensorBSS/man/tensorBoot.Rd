\name{tensorBoot}
\alias{tensorBoot}

\title{
Bootstrapping or Permuting a Data Tensor
}
\description{
The function takes bootstrap samples or permutes its content along the last dimension of the tensor.
}
\usage{
tensorBoot(x, replace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Array of an order of at least two with the last dimension corresponding to the sampling units.}
  \item{replace}{Logical. Should sampling be performed with or without replacement.
}
}
\details{
Assume an array of dimension \eqn{r+1}, where the last dimension represents the \eqn{n} sampling units and the first \eqn{r} dimensions the data per unit. The function then returns an array of the same dimension as \code{x} where either \eqn{n} bootstraps samples are selected or the units are permuted.
}
\value{
The bootstrapped or permuted samples in an array with the same dimension as \code{x}.
}

\author{
Christoph Koesner
}


\examples{
x <- array(1:50, c(2, 5, 5))
x
tensorBoot(x)
tensorBoot(x, replace = FALSE)

x <- array(1:100, c(2, 5, 2, 5))
x
tensorBoot(x)
}
\keyword{ utilities }
\keyword{ array }