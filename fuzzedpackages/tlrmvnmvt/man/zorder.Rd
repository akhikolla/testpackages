\name{zorder}
\alias{zorder}
\title{Index locations on the 2D plane}
\description{
  Index a set of locations of the 2D plane with a Z-curve, when 
  they are scaled into the unit square. The goal of this function 
  is to make locations close in geometry also close in the index. When
  partitioned into clusters, the inter-cluster correlation is more 
  likely to be low-rank.
}
\usage{
zorder(geom)
}
\arguments{
  \item{geom}{2D geometry, each row represents the location of a variable. 
    The geometry should be scaled into the unit square, (0,1) X (0,1)}
}
\value{
  A vector of indices based on the Z-curve.
}
\references{
  Foley, James, Andries van Dam, Steven Feiner, and John Hughes. 
    "Computer Graphics: Principle and Practice". Addison-Wesley 
    Publishing Company. Reading, Massachusetts: 1987. pages 870-871
}
\examples{
  n = 333
  set.seed(0)
  geom = matrix(runif(n*2), n, 2)
  idx = zorder(geom)
  idx
}
\author{Jian Cao, Marc Genton, David Keyes, George Turkiyyah}

