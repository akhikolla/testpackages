\name{tensorTransform2}
\alias{tensorTransform2}

\title{
Linear Transformations of Tensors from Several Modes
}

\description{
Applies a linear transformation to user selected modes of each individual tensor in an array of tensors. The function is a generalization of \code{\link{tensorTransform}} which only transforms one specific mode.
}

\usage{
tensorTransform2(x, A, mode, transpose = FALSE)
}

\arguments{
  \item{x}{Array of order r+1 >= 2 where the last dimension corresponds to the sampling units.}
  \item{A}{A list of r matrices to apply linearly to the corresponding mode.}
  \item{mode}{subsetting vector indicating which modes should be linearly transformed by multiplying them with the corresponding matrices from \code{A}.}
  \item{transpose}{logical. Should the matrices in \code{A} be transposed before the mode wise transformations or not.}
}
\details{
For the modes \eqn{i_1,\ldots,i_k}{i_1,...,i_k}, specified via \code{mode}, the function applies the linear transformation given by the matrix \eqn{A^{i_j}}{A^(i_j)} of size \eqn{q_{i_j} \times p_{i_j}}{q_i_j x p_i_j} to the \eqn{i_j}{i_j}th mode of each of the \eqn{n}{n} observed tensors \eqn{X_{i_j}}{X_i_j} in the given \eqn{p_1 \times p_2 \times \ldots \times p_r \times n}{p_1 x p_2 x ... x p_r x n}-dimensional array \code{x}. 
}

\value{
Array with r+1 dimensions where the dimensions specfied via \code{mode} are transformed.
}

\author{
Klaus Nordhausen
}

\seealso{
\code{\link{tensorTransform}}
}

\examples{
n <- 5
x <- array(rnorm(5*6*7), dim = c(7, 6, 5))
A1 <- matrix(runif(14), ncol = 7)
A2 <- matrix(rexp(18), ncol = 6)
A  <- list(A1 = A1, A2 = A2)
At <- list(tA1 = t(A1), tA2 = t(A2))

x1 <- tensorTransform2(x, A, 1)
x2 <- tensorTransform2(x, A, -2)
x3 <- tensorTransform(x, A1, 1)
x1 == x2
x1 == x3
x4 <- tensorTransform2(x,At,-2, TRUE)
x1 == x4
x5 <- tensorTransform2(x, A, 1:2)
}

\keyword{ utilities }
\keyword{ array }