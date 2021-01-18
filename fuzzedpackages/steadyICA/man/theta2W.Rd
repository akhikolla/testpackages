\name{theta2W}
\alias{theta2W}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert angles to an orthogonal matrix.
}
\description{
Convert d*(d-1)/2 angles from a sequence of Givens rotations to a d x d orthogonal matrix.
}
\usage{
theta2W(theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{A scalar or vector of length d*(d-1)/2 of values from which the d x d orthogonal matrix is calculated.
}
}

\value{
A d x d orthogonal matrix resulting from the sequence of d*(d-1)/2 Givens rotation matrices.
}
\author{
David S. Matteson}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
Golub, G. & Van Loan, C. 1996. Matrix computations. Johns Hopkins University Press.
}
\seealso{
\code{\link{W2theta}}
}
\examples{
#Generate orthogonal matrix:
mat <- matrix(rnorm(9),3,3)
W = svd(mat)$u

theta <- W2theta(W)

#Recovers W:
theta2W(theta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Givens }
\keyword{ orthogonal }% __ONLY ONE__ keyword per line
