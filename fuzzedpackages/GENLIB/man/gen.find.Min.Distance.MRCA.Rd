\name{gen.find.Min.Distance.MRCA}
\alias{gen.find.Min.Distance.MRCA}
\title{Finds the minimum distances between pairs of individuals given a set of ancestors.}
\description{Returns the shortest distances (number of meioses) between pairs of probands given the matrix of MRCAs output by the gen.findMRCA function.}
\usage{gen.find.Min.Distance.MRCA(genMatrix, individuals="ALL", ancestors="ALL")}
\arguments{
\item{genMatrix}{Matrix of most recent common ancestors, MRCAs, obtained with gen.findMRCA. Required.}
\item{individuals}{Vector of proband id numbers to include. All are included by default.}
\item{ancestors}{Vector of MRCA id numbers to include. All are included by default.}
}
\value{returns a matrix}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.founder}}
\code{\link{gen.findMRCA}}
\code{\link{gen.findDistance}}
\code{\link{gen.findFounders}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
matMRCA<-gen.findMRCA(genJi, individuals=c(1,29), NbProcess = 1)
gen.find.Min.Distance.MRCA(matMRCA)

\donttest{
# Increasing NbProcess will decrease execution time
data(genea140) 
gen140<-gen.genealogy(genea140) 
matMRCA<-gen.findMRCA(gen140, individuals=c(409033,408728,408828), NbProcess = 1)
gen.find.Min.Distance.MRCA(matMRCA, individuals=c(409033,408728))
}
}
\keyword{manip}
