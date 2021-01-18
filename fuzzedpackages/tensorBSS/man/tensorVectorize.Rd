\name{tensorVectorize}
\alias{tensorVectorize}

\title{Vectorize an Observation Tensor}
\description{
Vectorizes an array of array-valued observations into a matrix so that each column of the matrix corresponds to a single observational unit.
}
\usage{
tensorVectorize(x)
}
\arguments{
  \item{x}{Array of an order at least two with the last dimension corresponding to the sampling units.}
}
\details{
Vectorizes a \eqn{p_1 \times p_2 \times \ldots \times p_r \times n}{p_1 x p_2 x ... x p_r x n}-dimensional array into a \eqn{p_1 p_2 \ldots p_r \times n}{p_1 p_2 ... p_r x n}-dimensional matrix, each column of which then corresponds to a single observational unit. The vectorization is done so that the \eqn{r}{r}th index goes through its cycle the fastest and the first index the slowest.

Note that the output is a matrix of the size "number of variables" x "number of observations", that is, a transpose of the standard format for a data matrix. 
}
\value{
Matrix whose columns contain the vectorized observed tensors.
}
\author{
Joni Virta
}
\examples{
# Generate sample data.
n <- 100
x <- t(cbind(rnorm(n, mean = 0),
             rnorm(n, mean = 1),
             rnorm(n, mean = 2),
             rnorm(n, mean = 3),
             rnorm(n, mean = 4),
             rnorm(n, mean = 5)))
             
dim(x) <- c(3, 2, n)

# Matrix of vectorized observations.
vecx <- tensorVectorize(x)

# The covariance matrix of individual tensor elements
cov(t(vecx))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ array }
