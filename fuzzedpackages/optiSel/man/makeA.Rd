\name{makeA}
\Rdversion{1.1}
\alias{makeA}
\title{Calculates the Pedigree-based Additive Relationship Matrix 
}
\description{                                                 
Calculates the the Pedigree-based Additive Relationship Matrix. This is twice the pedigree based kinship matrix.
}
\usage{
makeA(Pedig, keep.only=NULL, keep=keep.only, AFounder=NULL)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree. The data frame has columns (1) Individual, (2) Sire, (3) Dam. Missing parents are coded as NA. Both parents must either be missing or present. If this is not the case use \link{prePed}.}
\item{keep}{If \code{keep} is provided then kinships are computed only for these animals and their ancestors.}
\item{keep.only}{If \code{keep.only} is provided then kinships are computed only for these animals.}
\item{AFounder}{Additive relationship matrix of the founders. The row names are the ids of the founders. By default, founders are assumed to be unrelated. Founders not included in this matrix are also assumed to be unrelated.}
}

\details{
Computation of pedigree based additive relationship matrix A which is twice the kinship matrix. For individuals i and j it is defined as
\tabular{l}{
 Aij = 2*(Probability that two alleles chosen from individuals i and j are IBD). \cr
}
}

\value{
Additive relationship matrix.
}


\examples{
data(PedigWithErrors)
data(Phen)
Pedig <- prePed(PedigWithErrors)
keep  <- Pedig$Indiv[summary(Pedig)$equiGen>5 & Pedig$Indiv \%in\% Phen$Indiv]
A     <- makeA(Pedig, keep.only=keep)
A[1:3,1:3]
}



\author{Robin Wellmann}
