\name{gen.findMRCA}
\alias{gen.findMRCA}
\title{Finding most recent common ancestors, MRCAs}
\description{Returns MRCAs of pairs of specified individuals and the distance (number of meioses) between individuals through the MRCAs.}
\usage{gen.findMRCA(gen, individuals, NbProcess=parallel::detectCores()-1)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{individuals}{Vector of individual id numbers for which to find the MRCAs. Required}
\item{NbProcess}{Number of processes to use when running this function. Default=parallel::detectCores()-1}
}
\value{returns a matrix}
\note{This function uses the parallel programming functions of foreach, snow and doSNOW for launching several processes.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.founder}}
\code{\link{gen.findFounders}}
\code{\link{gen.findDistance}}
\code{\link{gen.find.Min.Distance.MRCA}}
}
\examples{
data(geneaJi)
genJi<-gen.genealogy(geneaJi)
gen.findMRCA(genJi, individuals=c(1,29), NbProcess = 1)

\donttest{
 # For a more complex example:
 data(genea140) 
 gen140<-gen.genealogy(genea140) 
 all_commonFounders<-gen.findFounders(gen140, individuals=c(409033,408728), NbProcess = 1)
 # 127 founders common to #409033 and #408728
 length(all_commonFounders)
 # 22 most recent common ancestors for #409033 and #408728
 MRCA_2ind<-gen.findMRCA(gen140, individuals=c(409033,408728), NbProcess = 1)
}
}
\keyword{manip}
