\name{rtmg}
\alias{rtmg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample from truncated multivariate Gaussians
}
\description{
This function generates samples from a Markov chain whose equilibrium distribution is 
a d-dimensional multivariate Gaussian truncated by  linear and quadratic inequalities.
The probability log density is 

log p(X) = - 0.5 X^T M X + r^T X  + const

in terms of a precision matrix M and a vector r. The constraints are imposed  as explained below. The Markov chain is built using the Hamiltonian Monte Carlo technique. See the reference below for an explanation of the method.
}
\usage{
rtmg(n, M, r, initial, f = NULL, g = NULL, q = NULL, burn.in = 30)
}
\arguments{
  \item{n}{
 Number of samples. 
}
  \item{M}{
A d-by-d precision matrix of the multivariate Gaussian density. 
}
  \item{r}{
A d-dimensional vector for the linear coefficient of the log density.
}
  \item{initial}{
A d-dimensional vector with the initial value of the Markov chain. Must satisfy the truncation inequalities strictly.
}
  \item{f}{ An m-by-d matrix, where m is the number of linear constraints. The constraints require each component of the m-dimensional vector

fX + g

to be non-negative.
}
 \item{g}{ An m-dimensional vector with the constant terms in the above linear constraints.

}
  \item{q}{ A list of quadratic constraints. Each element i of q is a 3-element list of the form 

q[[i]] = list(A,B,C)
 
and imposes a quadratic constraint of the form

X^T A X + B^T X + C >=0

where A is a d-by-d matrix, B is a d-dimensional vector and C is a number.
}
  \item{burn.in}{
The number of burn-in iterations. The Markov chain is sampled n + burn.in times, and the last n samples are returned.
}
}
\details{
For linear constraints, both f and g must be provided. If no values for f, g and q are provided, the function returns samples from an untruncated Gaussian.
When the truncation equations define several disconnected regions, the samples belong to the region of the initial value.
}
\value{ An n-by-d matrix, where each row contains a sample.
}
\references{
Pakman, A. and Paninski, L., Exact Hamiltonian Monte Carlo for Truncated Multivariate Gaussians - Journal of Computational and Graphical Statistics, 2014

\url{http://arxiv.org/abs/1208.4118}
}
\author{
Ari Pakman

Maintainer: Ari Pakman <ari@stat.columbia.edu>
}
%\note{
%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# Set number of samples
n=15000;

#Define precision matrix and linear term 
M = matrix(c(.5,-.4, -.4,.5), 2,2)
r = c(0,0)

# Set initial point for the Markov chain
initial = c(4,1)

# Define two linear constraints
f = diag(2)
f[1,2] = 1
g = c(0,0)

# Define two quadratic constraints
A1 = matrix(c(-1/8,0,0,-1/2),2,2)
B1 = c(.5,.5)
C1 = 3/4
constr1 = list(A1,B1,C1)

A2 = matrix(c(4,-1,-1,8),2,2)
B2 = c(0,5)
C2 = -1
constr2 = list(A2,B2,C2)

q = list(constr1,constr2)

# Sample and plot
samples = rtmg(n, M, r, initial, f,g, q); 
plot(samples, pch=".")

}
\keyword{distribution } % __ONLY ONE__ keyword per line
