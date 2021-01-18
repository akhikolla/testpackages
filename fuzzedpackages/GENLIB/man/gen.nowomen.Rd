\name{gen.nowomen}
\alias{gen.nowomen}
\title{Number of women}
\description{Returns the number of women included in the genealogy}
\usage{gen.nowomen( gen)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
}
\value{returns a vector of integer}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.nomen}}
\code{\link{gen.nochildren}}
\code{\link{gen.noind}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gen.nowomen(genJi)

}
\keyword{utilities}
