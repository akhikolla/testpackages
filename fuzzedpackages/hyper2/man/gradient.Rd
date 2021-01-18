\name{gradient}
\alias{gradient}
\title{Differential calculus}
\description{
The gradient of the log-likelihood of a hyper2 object, at a specific
point in probability space
}
\usage{
gradient(H, probs)
}
\arguments{
  \item{H}{A \code{hyper2} object}
  \item{probs}{A vector of probabilities}
}
\details{

  Function \code{gradient()} returns the gradient of the log-likelihood
  function.  If the \code{hyper2} object is of size \eqn{n}, then
  argument \code{probs} must be a vector of length \eqn{n-1}.
}
\author{
Robin K. S. Hankin
}
\note{
This functionality is peculiarly susceptible to off-by-one errors.
}
\seealso{\code{\link{differentiate}}}
\examples{

data(chess)
p <- c(1/2,1/3)
delta <- rnorm(2)/1e5  # delta needs to be quite small

deltaL  <- loglik(chess,p+delta) - loglik(chess,p)
deltaLn <- sum(delta*gradient(chess,p + delta/2))   # numeric

deltaL - deltaLn  # should be small

}
