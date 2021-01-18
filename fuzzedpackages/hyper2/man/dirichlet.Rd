\name{dirichlet}
\alias{dirichlet}
\alias{gd}
\alias{Dirichlet}
\alias{gd}
\alias{GD}
\alias{GD_wong}
\alias{Connor}
\alias{Mosimann}
\title{
Dirichlet distribution and generalizations
}
\description{
The Dirichlet distribution in likelihood (for p) form, including the
generalized Dirichlet distribution due to Connor and Mosimann
}
\usage{
dirichlet(powers, alpha, pnames=NA)
GD(alpha, beta, beta0=0, pnames=NA)
GD_wong(alpha, beta, pnames=NA)
}
\arguments{
  \item{powers}{In function \code{dirichlet()} a (named) vector of powers}
  \item{alpha,beta}{A vector of parameters for the Dirichlet or
    generalized Dirichlet distribution}
  \item{beta0}{In function \code{GD()}, an arbitrary parameter}
  \item{pnames}{Optional; a character vector specifying the names of
    \eqn{p_1} through \eqn{p_n}}
}
\details{
These functions are really convenience functions.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{hyper2}}}
\examples{

x1 <- dirichlet(c(a=1,b=2,c=3))
x2 <- dirichlet(c(c=3,d=4))


# will not work: x1+x2  ## probably not what you wanted

}
