\name{whitener}
\alias{whitener}
\alias{whiten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Whitening function}
\description{
Subtract column means and transform columns such that the empirical covariance is equal to the identity matrix. Uses the SVD.}
\usage{
whitener(X, n.comp = ncol(X), center.row = FALSE, irlba = FALSE)
}
\arguments{
  \item{X}{
n x p matrix
}

  \item{n.comp}{
number of components to retain, i.e., first n.comp left eigenvectors from svd are retained
}

  \item{center.row}{
center both rows and columns prior to applying SVD (the resulting whitened data does not have zero-mean rows)
}
  \item{irlba}{if TRUE, uses irlba to approximate the first n.comp left eigenvectors. See Note.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{whitener}{the matrix such that X\%*\%whitener has zero mean and covariance equal to the identity matrix}
\item{Z}{the whitened data, i.e., X\%*\%whitener = Z}
}
\author{
Benjamin Risk}
\note{
The use of the option 'irlba = TRUE' requires the package irlba and is very useful
for large p. The function irlba only calculates
the first n.comp eigenvectors and is much faster than svd for p >> n.comp, for e.g., in groupICA of fMRI data.}
\seealso{
\code{\link{svd}}, 
\code{\link[irlba:irlba]{irlba::irlba}}
}
\examples{
simData <- cbind(rnorm(1000,1,2),rnorm(1000,-1,3),rnorm(1000,4,1))
simMVN <- simData\%*\%matrix(rnorm(12),3,4)
simWhiten <- whitener(simMVN,n.comp = 3)
colMeans(simWhiten$Z)
cov(simWhiten$Z)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ whiten}
\keyword{ prewhiten }
\keyword{ standardize }
