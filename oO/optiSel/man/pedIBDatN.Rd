\name{pedIBDatN}
\Rdversion{1.1}
\alias{pedIBDatN}
\title{Calculates the Pedigree Based Kinship at Native Alleles
}
\description{
Calculates the kinship at native alleles, which is the pedigree based probability of native alleles to be IBD.
}
\usage{
pedIBDatN(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only, nGen=NA)}

\arguments{
\item{Pedig}{
Data frame containing the pedigree with \code{Indiv} (Individual ID), \code{Sire}, and \code{Dam} in the first 3 columns, column \code{Breed} with breed names, and possibly column \code{Sex}. Missing parents are coded as NA, 0, or "0".
}
\item{thisBreed}{Name of the breed for which the kinships are to be computed.}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
\item{nGen}{Number of generations taken into account for estimating the native effective size. The default means that the native effective size is not estimated, which requires less memory.}
}

\details{
Calculates a list containing matrices needed to compute pedigree based kinships at native alleles, defined as the conditional probability that two randomly chosen alleles are IBD, given that both originate from native founders.
A native founder is an individual with unkown parents belonging to \code{thisBreed}.

The kinship at native alleles between individuals i and j is \code{Q1[i,j]/Q2[i,j]}.

The mean kinship at native alleles in the offspring is \code{(x'Q1x+d1)/(x'Q2x+d2)}, where \code{x} is the vector with genetic contributions of the selection candidates.

The native effective size is estimated from \code{nGen} generations only if \code{nGen} is not \code{NA}.
}

\value{
A list of class \code{ratioFun} including components:
\item{Q1}{matrix with \code{Q1[i,j]}  = Probability that two alleles chosen from individuals i and j are IBD and are native.}
\item{Q2}{matrix with \code{Q2[i,j]} = Probability that two alleles chosen from individuals i and j are both native.}
\item{d1}{The value by which the probability that two alleles chosen from the offspring are IBD and native increases due to genetic drift.}
\item{d2}{The value by which the probability that two alleles chosen from the offspring are native increases due to genetic drift.}
\item{id}{IDs of the individuals for which the probabilites have been computed.}
\item{mean}{Mean kinship at native alleles of the individuals specified in argument \code{keep.only}. Note that \code{1-mean} is the genetic diversity at native segments of the specified individuals from \code{thisBreed}}


}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
pKinatN <- pedIBDatN(Pedig, thisBreed="Hinterwaelder", keep.only=keep, nGen=6)

#Number of Migrant Founders: 237
#Number of Native  Founders: 150
#Individuals in Pedigree   : 1658
#Native effective size     : 49.5

## Mean kinship at native segments:
pKinatN$mean
#[1] 0.0776925

## Note that this can not be computed as mean(pKinatN$of).

## Results for individuals:
pKinatN$of <- pKinatN$Q1/pKinatN$Q2
pKinatN$of["276000812497583","276000812496823"]
#[1] 0.05941229


}



\author{Robin Wellmann}
