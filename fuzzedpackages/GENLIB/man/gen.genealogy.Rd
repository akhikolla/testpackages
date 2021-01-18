\name{gen.genealogy}
\alias{gen.genealogy}
\title{Create object of class GLgen}
\description{Creates an object of class GLgen that contains the ascending genealogies derived from input data in pedigree format}
\usage{gen.genealogy(ped, autoComplete=FALSE, ...)}
\arguments{
\item{ped}{A matrix or data frame with the following columns: individual id, father id, mother id, sex. Required.
	Individal id numbers must be numeric and unique.
	If an individual does not have a father and/or mother, the father and/or mother id numbers must be set to 0.
	Sex must be equal to M or 1 for males and F or 2 for females.
	The sex column is optional for this function but necessary for some other functions using GLgen objects.}
\item{autoComplete}{If TRUE, any non-zero father and mother id numbers not appearing in the individual id column, 
	will be added in the individual column as having no father or mother (both set to 0). Default to FALSE.}
\item{...}{Option to pass additionnal arguments automaticaly between methods. Internal use only.}
}
\value{An object of class GLgen.}
\seealso{
\code{\link{gen.graph}}
\code{\link{gen.genout}}
}
\examples{
ind<-c(1,2,3,11:21,101:108,201:202)
father<-c(11,15,15,102,0,102,0,103,103,103,105,105,107,107,0,202,0,202,202,0,202,0,0,0)
mother<-c(12,14,14,101,0,101,0,104,104,104,106,106,108,108,0,201,0,201,201,0,201,0,0,0)
sex<-c(1,2,2,1,2,1,2,1,2,1,2,2,2,1,2,1,1,2,1,2,1,2,2,1)
gen.df<-data.frame(ind, father, mother, sex)  
gen.genealogy(gen.df)

data(geneaJi) 
geneaJi[1:5,]
genJi<-gen.genealogy(geneaJi)
# Print basic information about the genealogy
genJi
}
\keyword{manip}
