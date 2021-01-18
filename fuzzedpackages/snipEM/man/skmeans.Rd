\name{skmeans}
\alias{skmeans}
\title{Snipped k-means clustering with cellwise outliers}
\description{
Perform k-means clustering on a data matrix with cellwise outliers using a snipping algorithm.
}
\usage{
skmeans(X, k, V, clust, s, itersmax = 10^5, D = 1e-1)
}
\arguments{
  \item{X}{Data.}
  \item{k}{Integer; number of clusters, \code{k>1}.}
  \item{V}{Binary matrix of the same size as X. Zeros correspond to initial snipped entries.}
  \item{clust}{Vector of size \code{n} containing values from \code{1} to \code{k}. Starting solution for class labels.}
  \item{itersmax}{Max number of iterations of the algorithm. Default is \code{3*10^5}.}
  \item{s}{Binary vector of size \code{n} for trimming, starting solution. Number of zeros will be preserved and correspond to trimmed rows. 
  If the vector is \code{rep(1,n)}, it performs no trimming. Default is \code{rep(1,n)}.}
  \item{D}{Tuning parameter for the fitting algorithm. Corresponds approximately to the maximal change in loss by switching two non
 outlying entries. Comparing different choices is recommended. Default is \code{1e-1}.}
}
\details{
  This function computes the \code{skmeans} estimator of Farcomeni
   (2014). It leads to robust k-means in presence of
   entry-wise and cellwise outliers. The number of snipped entries
   \code{sum(1-V)} and trimmed rows \code{sum(1-s)} is kept
   fixed throughout. Initial estimates for \code{V}, \code{s} and
   \code{clust} should be provided. Note that initializing with labels arising from
   classical (non-robust) clustering methods may be detrimental for the final
   performance of \code{skmeans} and may even yield an error due to
   empty clusters. 
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{loss} \tab Loss function (the total sum of squares) at convergence. \cr
	\code{mu} \tab Estimated locations. \cr
	\code{s} \tab Final (optimal) trimmed rows in vector of size \code{n}. \cr 
	\code{V} \tab Final (optimal) V matrix. \cr 
	\code{clust} \tab Final (optimal) class labels as vector of size \code{n}. \cr 
	}
}
\references{
Farcomeni, A. (2014) Snipping for robust k-means clustering under component-wise contamination, \emph{Statistics and Computing}, \bold{24}, 909-917
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\seealso{
  \code{\link{sclust}}, \code{\link{stEM}}, \code{\link{snipEM}},
  }
\examples{
set.seed(1234)
X <- matrix(NA,200,5)
# two clusters
k <- 2
X[1:100,] <- rnorm(100*5)
X[101:200,] <- rnorm(100*5,15)
clust <- rep(c(1,2), each=100)

# 5\% cellwise outliers
s <- sample(200*5,200*5*0.05)
X[s] <- runif(200*5*0.05,-100,100)
V <- X
V[s] <- 0
V[-s] <- 1

# Initial V and R
Vinit <- matrix(1, nrow(X), ncol(X))
Vinit[which(X > quantile(X,0.975) | X < quantile(X,0.025))] <- 0
km <- kmeans(X,k)
clustinit <- km$clust

# Snipped robust clustering
skm <- skmeans(X, k, Vinit, clustinit)

table(clust,km$clust)
table(clust,skm$clust)
}


