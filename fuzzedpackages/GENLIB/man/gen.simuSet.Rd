\name{gen.simuSet}
\alias{gen.simuSet}
\title{Gene dropping simulations with specified transmission probabilities}
\description{Returns the number of alleles transmitted to specified probands from ancestors considering sepcified transmission probabilities.}
\usage{gen.simuSet(gen, pro, ancestors, stateAncestors, 
                   probMatrix=matrix(c(
                                    c(1,0.5,0,0.5,0.25,0,0,0,0,1,1,1,1,0.75,0.5,1,0.5,0),
                                    c(1,0.5,0,0.5,0.25,0,0,0,0,1,1,1,1,0.75,0.5,1,0.5,0)),
                                     nrow=3, ncol=12),
                   simulNo = 5000)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Required.}
\item{ancestors}{Vector of ancestors id numbers to be included. Required.}
\item{stateAncestors}{Required vector indicating, for each ancestor in ancestors, the genotype state:

	0 = no disease allele present

	1 = 1 disease allele present

	2 = 2 disease alleles present}
\item{probMatrix}{Matrix of transmission probabilities of 0, 1, or 2 disease alleles. Default is:

	 [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]

[1,]  1.0 0.50    0    1 1.00  1.0  1.0 0.50    0     1  1.00   1.0

[2,]  0.5 0.25    0    1 0.75  0.5  0.5 0.25    0     1  0.75   0.5

[3,]  0.0 0.00    0    1 0.50  0.0  0.0 0.00    0     1  0.50   0.0

	where rows are genotype states (0,1,2) for the father and columns are states for the mother (in blocks of 3).

	1st block (col 1-3) are the probabilities of transmitting 0 allele to a son.

	2nd block (col 4-6) are the probabilities of transmitting 1 allele to a son.

	3rd block (col 7-9) are the probabilities of transmitting 0 allele to a daughter.

	4th block (col 10-12) are the probabilities of transmitting 1 allele to a daughter.
}
\item{simulNo}{Number of simulations to perform. Default is 5000.}
}
\value{A matrix with number of columns equal to the number of simulations and number of rows equal to the number of probands.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.simuProb}}
\code{\link{gen.simuSample}}
}
\examples{
# A case where only male subjects can receive alleles
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
onlyThroughMale<-matrix(c(c(1,0.5,0,0.5,0.25,0,0,0,0,1,1,1,1,0.75,0.5,1,0.5,0), rep(1,18)), 
                        nrow=3, ncol=12) 
gen.graph(genJi, indVarAffected=c(28,27,25,17,12,8,4,1), varAffected=c(28,27,25,17,12,8,4,1))
simu_1000a<-gen.simuSet(genJi, pro=c(1,28), ancestors=c(17,25), stateAncestors=c(1,1), 
                        simulNo = 1000,probMatrix=onlyThroughMale)
# Number of alleles received by probants
table(simu_1000a)
# Number of alleles received by all probants at each simulation
table(colSums(simu_1000a))
# Number of alleles received by each probant
table(simu_1000a[1,],simu_1000a[2,])

# A case where subjects are limited to one copy compared to what is normally expected
max1Allreceived<-matrix(rep(c(1,0.5,0,0.5,0.25,0,0,0,0,1,1,1,1,1,1,1,1,1), 2) , nrow=3, ncol=12)  
gen.graph(genJi)
simu_1000b<-gen.simuSet(genJi, pro=c(1,29), ancestors=25, stateAncestors=1, simulNo=10000,
                        probMatrix=max1Allreceived)
# Normal case matrix not changed
simu_1000original<-gen.simuSet(genJi, pro=c(1,29), ancestors=25, stateAncestors=1, simulNo=10000)
# Number of alleles received by each probant
table(simu_1000b) 
table(simu_1000original) 
# Number of alleles received by all probants at each simulation
table(colSums(simu_1000b)) 
table(colSums(simu_1000original)) 
# Number of alleles received by each probant
table(simu_1000b[1,],simu_1000b[2,])
table(simu_1000original[1,],simu_1000original[2,]) 
}
\keyword{datagen}
