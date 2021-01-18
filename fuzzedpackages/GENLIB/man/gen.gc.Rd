\name{gen.gc}
\alias{gen.gc}
\title{Genetic contribution of ancestors}
\description{Returns the genetic contribution of ancestors to the gene pool of sepcifiec probands}
\usage{gen.gc(gen, pro=0, ancestors=0, vctProb=c(0.5,0.5,0.5,0.5), typeCG="IND")}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Default is 0, which will select all individuals without children.}
\item{ancestors}{Vector of ancestors id numbers to be included. Default is 0, which will select all individuals without parents.}
\item{vctProb}{Vector of transmission probabilities. The first two values indicate the probabilities of transmission of a father to his son and daughter, respectively, and the following two values are the same for the mother. Default is vctProb=c(0.5,0.5,0.5,0.5).}
\item{typeCG}{IF typeCG="IND" (default), the genetic contribution from each ancestor is calculated for each proband.
	If typeCG="MEAN", the average (over all probands) genetic contributions of each ancestor is returned.
	If typeCG="PRODUCT", the product (over all probands) of genetic contributions is returned for each ancestor.
	If typeCG="TOTAL", the sum (over all probands) of genetic contributions is returned for each ancestor.
	If typeCG="CUMUL", ancestors are ranked in decreasing order of total contribution and cumulative contribution is returned.}
}
\value{A matrix with rows corresponding to probands and columns corrsponding to ancestors when typeCG="IND".
 	  For the other typeCG values, rows are ancestors and there is one column containing mean, product, total or cumulative values.}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.rec}}
\code{\link{gen.occ}}
\code{\link{gen.implex}}
\code{\link{gen.meangendepth}}
\code{\link{gen.completeness}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
gc<-gen.gc(genJi, pro=c(1,29), ancestors=c(3,6,10,12,14,16,18,20,26,28))
gc
gc_cum<-gen.gc(genJi, pro=c(1,29), ancestors=c(3,6,10,12,14,16,18,20,26,28), type="CUMUL")
gc_cum

data(genea140) 
gen140<-gen.genealogy(genea140) 
gc<-gen.gc(gen140, pro=c(454422,676521,677273,717537,717634,717709,868572), 
           ancestors=c(18305,18528,31114,18325))
gc
# Mother-daughter transmission only with probability=0.5
gc_MoLi<-gen.gc(gen140, pro=c(217891,302718,409282,802424,409682,443151),
		ancestors=c(18321,218231,296200,39066,18679,442607),vctProb=c(0,0,0,0.5))
gc_MoLi
}
\references{
Roberts DF. (1968) Genetic effects of population size reduction. Nature, 220, 1084-1088. 

O'Brien E, Jorde LB, Ronnlof B, Fellman JO, Eriksson AW. (1988) Founder effect and genetic disease in Sottunga, Finland. American Journal of Physical Anthropology, 77, 335-346.
}
\keyword{manip}

