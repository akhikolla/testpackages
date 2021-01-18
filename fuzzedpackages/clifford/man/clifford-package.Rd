\name{clifford-package}
\alias{clifford-package}
\docType{package}
\title{
\packageTitle{clifford}
}
\description{
\packageDescription{clifford}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{clifford}
\packageIndices{clifford}

}
\author{
\packageAuthor{clifford}

Maintainer: \packageMaintainer{clifford}
}
\references{
  \itemize{
    \item J. Snygg (2012).  \emph{A new approach to differential
      geometry using Clifford's geometric Algebra}, Birkhauser.  ISBN
    978-0-8176-8282-8
    \item D. Hestenes (1987). \emph{Clifford algebra to geometric
      calculus},  Kluwer. ISBN 90-277-1673-0
    \item C. Perwass (2009). \emph{Geometric algebra with applications
      in engineering}, Springer.  ISBN 978-3-540-89068-3
    \item D. Hildenbrand (2013). \emph{Foundations of geometric
      algebra computing}.  Springer, ISBN 978-3-642-31794-1
  }
}
\keyword{ package }
\seealso{\code{\link{clifford}}}
\examples{

as.1vector(1:4)

as.1vector(1:4) * rcliff()

# Following from Ablamowicz and Fauser (see vignette):
x <- clifford(list(1:3,c(1,5,7,8,10)),c(4,-10)) + 2
y <- clifford(list(c(1,2,3,7),c(1,5,6,8),c(1,4,6,7)),c(4,1,-3)) - 1
x*y # signature irrelevant



}
