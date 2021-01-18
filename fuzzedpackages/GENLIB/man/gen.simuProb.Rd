\name{gen.simuProb}
\alias{gen.simuProb}
\title{Gene dropping simulations - Probabilities}
\description{Returns the probabilities that specified probands inherit disease alleles from ancestors.}
\usage{gen.simuProb(gen, pro, statePro, ancestors, stateAncestors, simulNo=5000, 
                    probRecomb=c(0,0), probSurvival=1.0)}
\arguments{
\item{gen}{An object of class GLgen obtained with gen.genealogy, gen.lineages or gen.branching. Required.}
\item{pro}{Vector of proband id numbers to be included. Required.}
\item{statePro}{Required vector indicating, for each proband in pro, the probability to be calculated:

	0 = no disease allele is transmitted

	1 = 1 disease allele is transmitted

	2 = 2 disease alleles are transmitted

	3 = 1 or 2 disease alleles are transmitted}
\item{ancestors}{Vector of ancestors id numbers to be included. Required.}
\item{stateAncestors}{Required vector indicating, for each ancestor in ancestors, the genotype state:

	0 = no disease allele present

	1 = 1 disease allele present

	2 = 2 disease alleles present}
\item{simulNo}{Number of simulations to perform. Default is 5000.}
\item{probRecomb}{Recombination probabilities for males and females. Default is no recombination.}
\item{probSurvival}{Survival probability for homozygotes. Default is 1.}
}
\value{A list containing the following:

 	the joint probability of specified statePro for all probands,

 	the probability of specified statePro for each proband,

 	the probability that, 0, 1, ..., and all probands inherit the specified number of disease alleles.
}
\seealso{
\code{\link{gen.genealogy}}
\code{\link{gen.simuSample}}
\code{\link{gen.simuSet}}
\code{\link{gen.simuSampleFreq}}
}
\examples{
data(geneaJi) 
genJi<-gen.genealogy(geneaJi) 
# Probability that subjects 1 and 29 get 1 and 2 alleles from ancestors 20 and 25,
# that have themselves 2 and 1.
gen.simuProb(genJi, pro=c(1,29), statePro=c(1,2), ancestors=c(20,25), stateAncestors=c(2,1),
             simulNo=10000)
# Probability that subjects 1 and 29 get 1 segment from ancestors 25,
# knowing the segment has a male
# recombination rate of 0.02 and a female recombination rate of 0.04.
gen.simuProb(genJi, pro=c(1,29), statePro=c(1,1), ancestors=c(25), stateAncestors=c(1),
             simulNo=10000, probRecomb = c(0.02, 0.04))
# Probability that subjects 1 and 29 get 1 and 2 alleles from ancestors 20 and 25,
# that have themselves 2 and 1 and knowing that homozygous people have a survival rate
# of 0.50.
gen.simuProb(genJi, pro=c(1,29), statePro=c(1,2), ancestors=c(20,25), stateAncestors=c(2,1),
             simulNo=10000, probSurvival=0.5)

}
\keyword{datagen}
