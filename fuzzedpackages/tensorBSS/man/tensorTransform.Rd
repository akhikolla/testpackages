\name{tensorTransform}
\alias{tensorTransform}
\title{Linear Transformation of Tensors from mth Mode}
\description{
Applies a linear transformation to the mth mode of each individual tensor in an array of tensors.
}
\usage{
tensorTransform(x, A, m)
}

\arguments{
  \item{x}{Array of an order at least two with the last dimension corresponding to the sampling units.}
  \item{A}{Matrix corresponding to the desired linear transformation with the number of columns equal to the size of the \code{m}th dimension of \code{x}.}
  \item{m}{The mode from which the linear transform is to be applied.}
}
\details{
Applies the linear transformation given by the matrix \eqn{A}{A} of size \eqn{q_m \times p_m}{q_m x p_m} to the \eqn{m}{m}th mode of each of the \eqn{n}{n} observed tensors \eqn{X_i}{X_i} in the given \eqn{p_1 \times p_2 \times \ldots \times p_r \times n}{p_1 x p_2 x ... x p_r x n}-dimensional array \code{x}. This is equivalent to separately applying the linear transformation given by \eqn{A}{A} to each \eqn{m}{m}-mode vector of each \eqn{X_i}{X_i}. 
}
\value{
Array of size \eqn{p_1 \times p_2 \times \ldots \times p_{m-1} \times q_m \times p_{m+1} \times \ldots \times p_r \times n}{p_1 x p_2 x ... x p_(m-1) x q_m x p_(m+1) x ... x p_r x n}
}
\author{
Joni Virta
}
\examples{
# Generate sample data.
n <- 10
x <- t(cbind(rnorm(n, mean = 0),
             rnorm(n, mean = 1),
             rnorm(n, mean = 2),
             rnorm(n, mean = 3),
             rnorm(n, mean = 4),
             rnorm(n, mean = 5)))

dim(x) <- c(3, 2, n)

# Transform from the second mode
A <- matrix(c(2, 1, 0, 3), 2, 2)
z <- tensorTransform(x, A, 2)

# Compare
z[, , 1]
x[, , 1]\%*\%t(A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ array }
