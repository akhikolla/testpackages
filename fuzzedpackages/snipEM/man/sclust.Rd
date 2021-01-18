\name{sclust}
\alias{sclust}
\title{Snipping for robust model based clustering analysis with
cellwise outliers}
\description{
  Estimates a finite Gaussian mixture model optimized over a
  snipping set. 
}
\usage{
sclust(X, k, V, R, restr.fact=12, tol = 1e-04, maxiters = 100, 
          maxiters.S = 1000, print.it = FALSE) 
}
\arguments{
  \item{X}{Data.}
  \item{k}{Number of clusters}
  \item{V}{Binary matrix of the same size as \code{X}. Zeros correspond to
  initial snipped entries.}
  \item{R}{Initial guess for cluster labels, \code{1} to \code{k}.}
  \item{restr.fact}{Restriction factor, i.e., constraint on the condition number of all covariance matrices for each cluster. Default is 12.}
  \item{tol}{Tolerance for convergence. Default is \code{1e-4}.}
  \item{maxiters}{Maximum number of iterations for the SM algorithm. Default is \code{100}.}
  \item{maxiters.S}{Maximum number of iterations of the inner greedy snipping algorithm. Default is \code{1000}.}
  \item{print.it}{Logical; if TRUE, partial results are print. Default is \code{FALSE}.}
}
\details{
   This function computes the \code{sclust} estimator of Farcomeni
   (2014). It leads to robust mixture modeling in presence of entry-wise outliers. It is
   based on a classification-expectation-snip-maximize (CESM) algorithm. At the S step, the
   likelihood is optimized over the set of snipped entries, at the M
   step the location and scatter estimates are updated. The S step is
   based on a greedy algorithm, unlike the one proposed in Farcomeni
   (2014,2014a). The number of snipped entries \code{sum(1-V)} is kept
   fixed throughout. Note that initializing with labels arising from
   classical (non-robust) clustering methods may be detrimental for the final
   performance of \code{sclust} and may even yield an error due to
   empty clusters. 
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{R} \tab Final cluster labels. \cr
	\code{mu} \tab Estimated location matrix. \cr
	\code{S} \tab Array of estimated scatter matrices. \cr 
	\code{V} \tab Final (optimal) V matrix. \cr 
	\code{lik} \tab Gaussian log-likelihood at convergence. \cr 
	\code{iter} \tab Number of outer iterations before convergence. \cr 
	}
}
\references{
Farcomeni, A. (2014) Snipping for robust k-means clustering under component-wise contamination, \emph{Statistics and Computing}, \bold{24}, 909-917

Farcomeni, A. (2014) Robust constrained clustering in presence of entry-wise outliers, \emph{Technometrics}, \bold{56}, 102-111   
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\seealso{
  \code{\link{snipEM}}, \code{\link{stEM}},
    \code{\link{sumlog}},
  \code{\link{ldmvnorm}} 
}
\examples{
set.seed(1234)
X <- matrix(NA,200,5)
# two clusters
k <- 2
X[1:100,] <- rnorm(100*5)
X[101:200,] <- rnorm(100*5,15)
R <- rep(c(1,2), each=100)

# 5\% cellwise outliers
s <- sample(200*5,200*5*0.05)
X[s] <- runif(200*5*0.05,-100,100)
V <- X
V[s] <- 0
V[-s] <- 1

# Initial V and R
Vinit <- matrix(1, nrow(X), ncol(X))
Vinit[which(X > quantile(X,0.975) | X < quantile(X,0.025))] <- 0
Rinit <- kmeans(X,2)$clust

# Snipped robust clustering
sc <- sclust(X,2,Vinit,Rinit)
table(R,Rinit)
table(R,sc$R)
}


