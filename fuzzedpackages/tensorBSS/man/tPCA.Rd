\name{tPCA}
\alias{tPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PCA for Tensor-Valued Observations}
\description{
Computes the tensorial principal components. 
}
\usage{
tPCA(x, p = NULL, d = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric array of an order at least three. It is assumed that the last dimension corresponds to the sampling units.}
  \item{p}{A vector of the percentages of variation per each mode the principal components should explain.}
  \item{d}{A vector of the exact number of components retained per each mode. At most one of this and the previous argument should be supplied.} 
}
\details{
The observed tensors (array) \eqn{X} of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} measured on \eqn{N} units are projected from each mode on the eigenspaces of the \eqn{m}-mode covariance matrices of the corresponding modes. As in regular PCA, by retaining only some subsets of these projections (indices) with respective sizes \eqn{d_1, d_2, ... d_r}, a dimension reduction can be carried out, resulting into observations tensors of size \eqn{d_1 \times d_2 \times \ldots \times d_r}{d_1 x d_2 x ... x d_r}. In R the sample of \eqn{X} is saved as an \code{\link{array}} of dimensions
\eqn{p_1, p_2, \ldots, p_r, N}{p_1, p_2, ..., p_r, N}.

}
\value{
A list containing the following components: 
\item{S }{Array of the same size as x containing the principal components.}
\item{U }{List containing the rotation matrices}
\item{D }{List containing the amounts of variance explained by each index in each mode.}
\item{p_comp}{The percentages of variation per each mode that the principal components explain.}
\item{Xmu }{The data location.}
}
\references{
\cite{Virta, J., Taskinen, S. and Nordhausen, K. (2016), Applying fully tensorial ICA to fMRI data, Signal Processing in Medicine and Biology Symposium (SPMB), 2016 IEEE,\cr \doi{10.1109/SPMB.2016.7846858}} 
}
\author{
Joni Virta
}


\examples{
# Digit data example

data(zip.train)
x <- zip.train

rows <- which(x[, 1] == 0 | x[, 1] == 1)
x0 <- x[rows, 2:257]
y0 <- x[rows, 1] + 1

x0 <- t(x0)
dim(x0) <- c(16, 16, 2199)


tpca <- tPCA(x0, d = c(2, 2))
pairs(t(apply(tpca$S, 3, c)), col=y0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{array }
