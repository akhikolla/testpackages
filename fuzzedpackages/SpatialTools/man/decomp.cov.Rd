\name{decomp.cov} 

\alias{decomp.cov}

\title{
	Calculates decomposition of covariance matrix
} 

\description{
	Calculates a decomposition of the provided covariance matrix, \code{V}, using the chosen method.  
} 

\usage{
	decomp.cov(V, method = "eigen")
} 

\arguments{
	\item{V}{A (symmetric, positive-definite) covariance matrix.} 
	\item{method}{A character vector specifying the method used to decompose \code{V}.  Options are "eigen", "chol", or "svd" (Eigen decomposition, Cholesky decomposition, or Singular value decomposition, respectively).} 
} 

\details{ The matrix \code{V} is assumed to be symmetric and positive definite.  Symmetry is checked, but the positive definiteness of the matrix is not.  Returns a decomposition matrix \code{U} such that \code{V} = \code{U} \%*\% \code{t(U)}.
}

\value{
	Returns a decomposition matrix \code{U} such that \code{V} = \code{U} \%*\% \code{t(U)}. 
}


%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{covariance} 
\keyword{decomposition} 

\seealso{
 \code{cov.sp}
} 

\examples{
	data(toydata)
	U <- decomp.cov(toydata$V, method = "chol")
	#range(toydata$V - U %*% t(U))
}