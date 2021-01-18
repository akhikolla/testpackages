\name{gen.genout}
\alias{gen.genout}
\title{Create pedigree data}
\description{Creates data frame in pedigree format from an object of class GLgen}
\usage{gen.genout( gen, sorted = "F")}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{sorted}{If sorted="F" (default), individual id number are not sorted in output data. Id numbers for parents will be placed before their children ids.
	If sorted="T", individual id numbers will be sorted.}
}
\value{returns a data.frame containing the following:
\strong{ind} 
\strong{father} 
\strong{mother} 
\strong{sex} 
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.graph}}
\code{\link{gen.branching}}
\code{\link{gen.lineages}}
}
\examples{
data(geneaJi) 
# Original data is a data.frame
geneaJi[1:12,]

genJi<-gen.genealogy(geneaJi) 
# As a genealogy object
genJi

# Genealogy as a data.frame
genJi_df<-gen.genout(genJi)
genJi_df[1:12,]

# Maternal lineage
genJi_MaLi<-gen.lineages(geneaJi, maternal = TRUE) 
# Maternal lineage as a data.frame
genJi_MaLi_df<-gen.genout(genJi_MaLi)
genJi_MaLi_df
}
\keyword{manip}
