\name{gen.noind}
\alias{gen.noind}
\title{Number of individuals}
\description{Returns the number of individuals included in the genealogy}
\usage{gen.noind( gen)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
}
\value{returns a vector of integer}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.nowomen}}
\code{\link{gen.nochildren}}
\code{\link{gen.nomen}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.noind(genJi)
}
\keyword{utilities}
