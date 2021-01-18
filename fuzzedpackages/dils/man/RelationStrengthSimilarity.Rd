\name{RelationStrengthSimilarity}
\alias{RelationStrengthSimilarity}
\title{Calculate the RSS from one node to another.}
\usage{
  RelationStrengthSimilarity(xadj, v1, v2, radius = 3,
    directed = TRUE,
    method = c("Rcpp", "BetterR", "NaiveR"))
}
\arguments{
  \item{xadj}{numeric matrix, then description of
  \code{arg1}.}

  \item{v1}{numeric Object type, then description of
  \code{arg2}.}

  \item{v2}{numeric Object type, then description of
  \code{arg2}.}

  \item{radius}{numeric, length of longest path examined
  from \code{v1} to \code{v2}.}

  \item{directed}{logical, if TRUE returns a symmetric RSS
  matrix.}

  \item{method}{character, choose the method of
  calculation.}
}
\value{
  numeric, Relation Strength Similarity score(s).
}
\description{
  For a single pair of nodes, implement the RSS algorithm
  of Chen et al. (2012).
}
\details{
  If \code{v1} and \code{v2} are specified, this returns
  the RSS from \code{v1} to \code{v2}.  If not, it
  calculates the RSS scores for all dyads in the network.
}
\examples{
g1 <- graph.atlas(128)
\dontrun{plot(g1)}

M1 <- as.matrix(get.adjacency(g1))
M1
RelationStrengthSimilarity(xadj=M1, v1=5, v2=6, radius=1)
RelationStrengthSimilarity(xadj=M1, v1=5, v2=6, radius=2)
RelationStrengthSimilarity(xadj=M1, v1=5, v2=6, radius=3)
RelationStrengthSimilarity(xadj=M1, v1=5, v2=6, radius=4)

RelationStrengthSimilarity(xadj=M1, radius=2)

TestUndirectedNetwork <- function(n) {
  M <- matrix(runif(n*n), nrow=n)
  M <- (M + t(M)) / 2
  diag(M) <- 0
  return(M)
}
M2 <- TestUndirectedNetwork(75)
system.time(RelationStrengthSimilarity(xadj=M2, directed=FALSE, method="BetterR"))  # all R
system.time(RelationStrengthSimilarity(xadj=M2, directed=FALSE))                    # Rcpp
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  "Discovering Missing Links in Networks Using Similarity
  Measures", Hung-Hsuan Chen, Liang Gou, Xiaolong (Luke)
  Zhang, C. Lee Giles. 2012.

  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{ScalablePCA}}
}

