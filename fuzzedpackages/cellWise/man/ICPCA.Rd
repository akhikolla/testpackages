\name{ICPCA}
\alias{ICPCA}

\title{
Iterative Classical PCA
}
\description{
This function carries out classical PCA when the data may contain
missing values, by an iterative algorithm. It is based on a Matlab function from the Missing Data Imputation Toolbox v1.0
 by A. Folch-Fortuny, F. Arteaga and A. Ferrer.
}
\usage{
ICPCA(X, k, scale = FALSE, maxiter = 20, tol = 0.005,
                 tolProb = 0.99, distprob = 0.99) 
}

\arguments{
  \item{X}{the input data, which must be a matrix or a data frame.
      It may contain NA's. It must always be provided.}
  \item{k}{the desired number of principal components}
  \item{scale}{a value indicating whether and how the original
              variables should be scaled. If \code{scale=FALSE} (default)
              or \code{scale=NULL} no scaling is performed (and a vector
              of 1s is returned in the \code{$scaleX} slot).
              If \code{scale=TRUE} the variables are scaled to have a
              standard deviation of 1. Alternatively scale can be a function like mad,
              or a vector of length equal to the number of columns
              of x. The resulting scale estimates are returned in the
              \code{$scaleX} slot of the output.}
  \item{maxiter}{ maximum number of iterations. Default is 20.}
  \item{tol}{ tolerance for iterations. Default is 0.005.}
  \item{tolProb}{tolerance probability for residuals. Defaults to 0.99.}
  \item{distprob}{probability determining the cutoff values for
              orthogonal and score distances. Default is 0.99.}
}

\value{
 A list with components: \cr
 \item{scaleX}{the scales of the columns of X.}
  \item{k}{the number of principal components.}
  \item{loadings}{the columns are the k loading vectors.}
  \item{eigenvalues}{the k eigenvalues.}
  \item{center}{vector with the fitted center.}
 \item{covmatrix}{estimated covariance matrix.}
  \item{It}{number of iteration steps.}
  \item{diff}{convergence criterion.}
  \item{X.NAimp}{data with all NA's imputed.}
  \item{scores}{scores of X.NAimp.}
  \item{OD}{orthogonal distances of the rows of X.NAimp.}
  \item{cutoffOD}{cutoff value for the OD.}
  \item{SD}{score distances of the rows of X.NAimp.}
  \item{cutoffSD}{cutoff value for the SD.}
  \item{indrows}{row numbers of rowwise outliers.}
  \item{residScale}{scale of the residuals.}
  \item{stdResid}{standardized residuals. Note that these are NA
                for all missing values of \code{X}.}
  \item{indcells}{indices of cellwise outliers.}
}
\references{
Folch-Fortuny, A., Arteaga, F., Ferrer, A. (2016). Missing Data Imputation Toolbox for MATLAB.
\emph{Chemometrics and Intelligent Laboratory Systems}, \bold{154}, 93-100.
}
\author{
Wannes Van Den Bossche
}

\examples{
library(MASS) 
set.seed(12345) 
n <- 100; d <- 10
A <- diag(d) * 0.1 + 0.9
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 100, FALSE)] <- NA
ICPCA.out <- ICPCA(x, k = 2)
plot(ICPCA.out$scores)
}
