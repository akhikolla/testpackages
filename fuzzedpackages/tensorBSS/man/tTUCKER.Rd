\name{tTUCKER}
\alias{tTUCKER}

\title{
Tucker (2) Transformation for a Tensor
}
\description{
This is a Tucker (2) transformation of a data tensor where the sampling dimension is uncompressed. The transfromation is known also under many different names like multilinear principal components analysis or generalized low rank approximation of matrices if the tensorial data is matrixvalued.
}
\usage{
tTUCKER(x, ranks, maxiter = 1000, eps = 1e-06)
}

\arguments{
  \item{x}{array with \eqn{r+1} dimensions where the last dimension corresponds to the sampling units.}
  \item{ranks}{vector of length r giving the dimensions of the compressed core tensor.}
  \item{maxiter}{maximum number of iterations for the algorithm.}
  \item{eps}{convergence tolerance.}
}
\details{
As initial solution \code{\link{tPCA}} is used and iterated using an alternating least squares (ALS) approach, known also as higher order orthogonal iteration (HOOI).
}

\value{
A list containing the following components: 
\item{S }{array of the compressed tensor.}
\item{U }{list containing the rotation matrices.}
\item{Xmu  }{the data location.}
\item{norm2xc}{squared norm of the original data tensor after centering.}
\item{norm2rxc}{squared norm of the reconstructed (centered) data tensor.}
\item{norm2ratio}{the ratio norm2rxc/norm2xc.}
\item{mEV}{list containing the eigenvalues from the m-mode covariance matrix when all but the relevant mode have be compressed.}
\item{tPCA }{The output from \code{\link{tPCA}} which was used as initial value.}
}
\references{
\cite{Lu, H., Plataniotis, K. and  Venetsanopoulos, A. (2008), MPCA: Multilinear principal component analysis of tensor objects, IEEE Transactions on Neural Networks, 19, 18-39. \doi{10.1109/TNN.2007.901277}}

\cite{Lietzen, N., Nordhausen, K. and Virta, J. (2019), Statistical analysis of second-order tensor decompositions, manuscript.} 
}

\author{
Klaus Nordhausen
}

\seealso{
\code{\link{tPCA}}
}

\examples{
data(zip.train)
x <- zip.train

rows <- which(x[, 1] == 0 | x[, 1] == 1)
x0 <- x[rows, 2:257]
y0 <- x[rows, 1] + 1

x0 <- t(x0)
dim(x0) <- c(16, 16, 2199)

tucker <-  tTUCKER(x0, ranks = c(2, 2), eps=1e-03)
pairs(t(apply(tucker$S, 3, c)), col=y0)

# To approximate the original data one uses then
x0r <- tensorTransform2(tucker$S, tucker$U)


}

\keyword{multivariate}
\keyword{array }