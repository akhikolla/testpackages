\name{pedIBD}
\Rdversion{1.1}
\alias{pedIBD}
\title{Calculates the Pedigree-based Kinship Matrix 
}
\description{                                                 
Calculates the \bold{ped}igree based probability of alleles to be \bold{IBD}. This pedigree based kinship matrix is also called coancestry matrix and is half the additive relationship matrix.
}
\usage{
pedIBD(Pedig, keep.only=NULL, keep=keep.only, kinFounder=NULL)}

\arguments{
\item{Pedig}{Data frame containing the pedigree with \code{Indiv} (individual ID), \code{Sire}, and \code{Dam} in the first 3 columns.
Missing parents are coded as NA. Both parents must either be missing or present. If this is not the case use function \link{prePed} to prepare the pedigree.}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
\item{kinFounder}{Kinship matrix for the founders. The row names are the ids of the founders. By default, founders are assumed to be unrelated. Founders not included in this matrix are also assumed to be unrelated.}
}

\details{
Computation of pedigree based kinship matrix f which is half the additive relationship matrix. For individuals i and j it is defined as
\tabular{l}{
 fij = Probability that two alleles chosen from individuals i and j are IBD. \cr
}
}

\value{
Kinship matrix.
}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv
Pedig <- prePed(PedigWithErrors, keep=keep, thisBreed="Hinterwaelder", lastNative=1970)
pedA  <- pedIBD(Pedig, keep.only=keep)
}



\author{Robin Wellmann}
