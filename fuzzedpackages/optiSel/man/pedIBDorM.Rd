\name{pedIBDorM}
\Rdversion{1.1}
\alias{pedIBDorM}
\title{Calculates Kinships taking Allele Origin into Account
}
\description{
Calculates the \bold{ped}igree based probability of alleles to be \bold{IBD} (identical by descent) \bold{or} \bold{M}igrant alleles: 
 For each pair of individuals the probability is computed that two alleles taken at random are IBD or are migrant alleles.
}
\usage{
pedIBDorM(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam, (4) Sex, and (5) Breed. Missing parents are coded as NA. Both parents must either be missing or present. If this is not the case use \link{prePed}.}
\item{thisBreed}{Name of the breed in column (5) of the pedigree for which the kinships are to be computed.}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
}
 
\details{
Computation of modified pedigree based kinship matrices taking allele origin into account. 

A native founder is an individual with unkown parents belonging to \code{thisBreed}. A migrant is an individual with unkown parents not belonging to \code{thisBreed}. 
}

\value{
A list with the following components:

\item{pedIBDorM}{Matrix containing for  individuals i and j the probability that two alleles chosen from the individuals are IBD or at least one of them is a migrant allele (only computed if 1 is in \code{method})}

\item{pedIBDorMM}{Matrix containing for  individuals i and j the probability that two alleles chosen from the individuals are IBD or both are migrant alleles (only computed if 2 is in \code{method})}

}

\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
Kin   <- pedIBDorM(Pedig, thisBreed="Hinterwaelder", keep.only=keep)

mean(Kin$pedIBDorM)
#[1] 0.8201792
mean(Kin$pedIBDorMM)
#[1] 0.335358
}

\author{Robin Wellmann}
