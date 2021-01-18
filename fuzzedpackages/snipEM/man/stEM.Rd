\name{stEM}
\alias{stEM}
\title{Snipping and trimming for location and scatter estimation with
casewise and cellwise outliers}
\description{
  Computes an estimator optimizing the Gaussian likelihood over a
  snipping and trimming set. 
}
\usage{
stEM(X, V, tol = 1e-4, maxiters = 500, maxiters.S = 1000, print.it = FALSE) 
}
\arguments{
  \item{X}{Data.}
  \item{V}{Binary matrix of the same size as \code{X}. Zeros correspond to
  initial snipped entries, rows of zeros correspond to initial trimmed
  entries.}
  \item{tol}{Tolerance for convergence. Default is \code{1e-4}.}
  \item{maxiters}{Maximum number of iterations for the SM algorithm. Default is \code{500}.}
  \item{maxiters.S}{Maximum number of iterations of the inner greedy snipping algorithm. Default is \code{1000}.}
   \item{print.it}{Logical; if \code{TRUE}, partial results are print. Default is \code{FALSE}.}
}
\details{
   This function combines computes the \code{snipEM} estimator of Farcomeni
   (2014) with trimming. Optimization over a trimming set is performed
   via usual concentration steps (Rousseeuw and van Driessen, 1999). 
   It therefore provides a robust estimate of
   location and scatter in presence of entry-wise and case-wise
   outliers. The number of snipped entries and trimmed rows is kept
   fixed throughout. \code{V} must contain at least one row of zeros
   (otherwise use \code{\link{snipEM}}). 
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{mu} \tab Estimated location. \cr
	\code{S} \tab Estimated scatter matrix. \cr 
	\code{V} \tab Final (optimal) V matrix. \cr 
	\code{lik} \tab Gaussian log-likelihood at convergence. \cr 
	\code{iter} \tab Number of outer iterations before convergence. \cr 
	}
}
\references{
Farcomeni, A. (2014) Snipping for robust k-means clustering under component-wise contamination, \emph{Statistics and Computing}, \bold{24}, 909-917

Farcomeni, A. (2014) Robust constrained clustering in presence of entry-wise outliers, \emph{Technometrics}, \bold{56}, 102-111   

Rousseeuw, P. J. and Van Driessen, K. (1999) A fast algorithm for the
minimum covariance determinant estimator, \emph{Technometrics}, \bold{41}, 212-223.
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\seealso{
  \code{\link{sclust}}, \code{\link{snipEM}},
    \code{\link{sumlog}},
  \code{\link{ldmvnorm}} 
}
\examples{
set.seed(1234)
X=matrix(rnorm(100*10),100,5)
X[1:5,]=50
X[6,1]=150

# initial V
V <- matrix(1, 100, 5)
V[1:5,]=0
Vtmp <- V[-c(1:5),]

# identify cells to be snipped
Vtmp[!is.na(match(as.vector(X[-c(1:5),]),boxplot(as.vector(X[-c(1:5),]),plot=FALSE)$out))] <- 0
V[-c(1:5),] <- Vtmp

resSTEM <- stEM(X, V)

}
