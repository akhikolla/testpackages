\name{noffspring}
\Rdversion{1.1}
\alias{noffspring}
\title{Calculates Optimum Numbers of Offspring}

\description{
Calculates the optimum numbers of offspring from optimum contributions of selection candidates.
}
\usage{
noffspring(cand, N, random=TRUE)}

\arguments{
\item{cand}{Data frame with optimum contributions (column \code{oc}), sexes (column \code{Sex}), and IDs (column \code{Indiv}) of the selection candidates.}
\item{N}{Desired number of individuals in the offspring population.}
\item{random}{Logical. If \code{2*N*oc[i]} is not an integer value (say \code{2*N*oc[i]=11.4}) then individual \code{i} will have either 11 or 12 offspring. The actual number is either determined randomly or not.}
}

\details{
The function calculates the optimum numbers of offspring of the selection candidates from the optimum contributions \code{cand$oc} and the size \code{N} of the offspring population.
}

\value{

Data frame with column \code{Indiv} containing the individual IDs and column \code{nOff} containing the optimum numbers of offspring. 

Column \code{nOff} is approximately \code{2*N*cand$oc} with \code{sum(noff[cand$Sex=="male"])=N} and \code{sum(noff[cand$Sex=="female"])=N}.
}


\examples{
set.seed(1)
data(PedigWithErrors)

Pedig      <- prePed(PedigWithErrors, thisBreed="Hinterwaelder")
use        <- Pedig$Born \%in\% (1998:2008) & Pedig$Breed=="Hinterwaelder"
Population <- sampleIndiv(Pedig[use, ], each=50)
pKin       <- pedIBD(Pedig, keep.only=Population)
Phen       <- Pedig[Population, ]
Phen$isCandidate <- Phen$Born \%in\% (2003:2008)

cont       <- agecont(Pedig, Population)
cand       <- candes(phen=Phen, fA=pedIBD(Pedig, keep.only=Phen$Indiv), cont=cont)
con        <- list(ub.fA=0.0175, uniform="female")
Offspring  <- opticont("max.BV", cand, con, trace = FALSE)

N <- 250
Candidate <- Offspring$parent
Candidate$nOff <- noffspring(Candidate, N)$nOff

sum(Candidate$nOff[Candidate$Sex=="male"])
#[1] 250

sum(Candidate$nOff[Candidate$Sex=="female"])
#[1] 250

round(2*N*Candidate$oc-Candidate$nOff, 2)
}

\author{Robin Wellmann}
