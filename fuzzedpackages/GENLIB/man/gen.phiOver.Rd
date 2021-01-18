\name{gen.phiOver}
\alias{gen.phiOver}
\title{Kinship above threshold}
\description{Returns the pairs of individuals with kinship coefficient values greater than specified threshold.}
\usage{gen.phiOver( phiMatrix, threshold)}
\arguments{
\item{phiMatrix}{A square matrix of kinship coefficients obtained for example with gen.phi. Required.}
\item{threshold}{Threshold of kinship values to return.}
}
\value{A data frame containing the probands and their kinship.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.phi}}
\code{\link{gen.phiMean}}
\code{\link{gen.phiCI}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
kinship<-gen.phi(genJi)
gen.phiOver(kinship, 0.1)

\donttest{
data(genea140) 
gen140<-gen.genealogy(genea140) 
phi9subj<-gen.phi(gen140,pro=c(408758,408950,409082,409111,676521,717537,717634,717709,868572))
gen.phiOver(phi9subj,0.025)
}
}
\keyword{manip}
