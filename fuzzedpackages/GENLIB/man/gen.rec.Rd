\name{gen.rec}
\alias{gen.rec}
\title{Ancestors coverage}
\description{Returns the number of specified probands genealogically related to specified ancestors (i.e., ancestor occurs in the proband's genealogy).}
\usage{gen.rec( gen, pro = "0", ancestors = "0", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{ancestors}{Vector of ancestors id numbers to be included. Default is 0, which will select all individuals without parents.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{A matrix with number of lines equal to the number of ancestors specified and one column.}
\note{If an ancestor is also a proband, he/she will be counted in his/her coverage.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.occ}}
\code{\link{gen.implex}}
\code{\link{gen.meangendepth}}
\code{\link{gen.gc}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
# Number of probands which are descendants of an ancestor
gen.rec(genJi)


data(genea140) 
gen140<-gen.genealogy(genea140) 
gen.rec(gen140, ancestors=c(18305, 18528, 31114, 18325))
}
\keyword{manip}
