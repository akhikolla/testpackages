\name{W2theta}
\alias{W2theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert an orthogonal matrix to its angular parameterization. }
\description{
Convert a d x d orthogonal matrix to a sequence of d*(d-1)/2 Givens rotations.
}
\usage{
W2theta(W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
  A d x d orthogonal matrix.
}
}
\details{
A d x d orthogonal matrix can be decomposed into a series of d*(d-1)/2 Givens rotation matrices, where each matrix is parameterized by a single angle.
}
\value{
A vector of length d*(d-1)/2 comprised of the angles.
}
\author{
David S. Matteson
}

\references{
Golub, G. & Van Loan, C. 1996. Matrix computations. Johns Hopkins University Press.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{theta2W}}
}
\examples{
theta = c(pi/6,pi/4,pi/2)

(W = theta2W(theta))

#Recover theta:
W2theta(W)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Givens }
\keyword{ orthogonal }% __ONLY ONE__ keyword per line
