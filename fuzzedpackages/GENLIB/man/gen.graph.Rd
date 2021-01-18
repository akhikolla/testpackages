\name{gen.graph}
\alias{gen.graph}
\title{Pedigree graphical tool}
\description{Function that plots pedigrees of GLgen objects.}
\usage{gen.graph( gen, pro=gen.pro(gen), ancestors=gen.founder(gen), 
       indVarAffected=gen.genout(gen)$ind, varAffected=gen.genout(gen)$ind, cex="1",
       col="0", symbolsize="1", width="1", packed="F", align="T", ...)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 'gen.pro(gen)', which will select all individuals without children.}
\item{ancestors}{Vector of ancestors id numbers to be included. Default is 'gen.founder(gen)', which will select all individuals without parents.}
\item{indVarAffected}{Vector of individuals id numbers used for labeling. Default is 'gen.genout(gen)$ind', which is all the individuals of the genealogy.}
\item{varAffected}{Vector of individuals labels. Default is 'gen.genout(gen)$ind', meaning that the individuals labels are their ids.}
\item{cex}{Controls text size (same as kinship2::plot.pedigree).  Default is 1.}
\item{col}{Color for each id (same as kinship2::plot.pedigree). Default is 0, which assigns the same color to everyone.}
\item{symbolsize}{Controls symbolsize (same as kinship2::plot.pedigree). Default is 1.}
\item{width}{For a packed pedigree, the minimum width allowed in the realignment of pedigrees (same as kinship2::plot.pedigree). Default is 1.}
\item{packed}{If T, uniform distance between all individuals at a given level (same as kinship2::plot.pedigree). Default is FALSE.}
\item{align}{Controls the extra effort spent trying to align children underneath parents, but without making the
          pedigree too wide (same as kinship2::plot.pedigree). Default is TRUE.}
\item{...}{Other arguments that can be passed to the kinship2::plot.pedigree.}
}
\value{Returns the same invisible list as that returned by kinship2::plot.pedigree, which contains the following:
\strong{plist} list containing the information about the pedigree (n, id, pos, fam, spouse)
\strong{x} x-axis position
\strong{y} y-axis position
\strong{boxw} box width
\strong{boxh} box height
\strong{call} the call made to plot.pedigree() (kinship2 package)
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.genout}}
\code{\link{gen.branching}}
}
\examples{
ind<-c(1,2,3,11:21,101:108,201:202)
father<-c(11,15,15,102,0,102,0,103,103,0,105,105,107,107,0,202,0,202,202,0,202,0,0,0)
mother<-c(12,14,14,101,0,101,0,104,104,0,106,106,108,108,0,201,0,201,201,0,201,0,0,0)
sex<-c(1,2,2,1,2,1,2,1,2,1,2,2,2,1,2,1,1,2,1,2,1,2,2,1)
gen.df<-data.frame(ind, father, mother, sex)  
genEx<-gen.genealogy(gen.df)
gen.graph(genEx)

data(geneaJi) 
geneaJi[1:5,]
genJi<-gen.genealogy(geneaJi)
gen.graph(genJi)
}
\keyword{aplot}
