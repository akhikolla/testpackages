\name{pogc}
\alias{pogc}
\title{
Parent Offspring Group Constructor
}
\description{
Assign offsprings to the parents.
}
\usage{
pogc(oh, genotypeError)
}
\arguments{
  \item{oh}{ \code{integer} opposing homozygotes matrix (Output of \code{\link{ohg}})
}
  \item{genotypeError}{\code{integer}
number of genotypeing error allowed in the \code{oh} matrix
}
}

\value{
Return a data frame with two columns. The first column is the animal ID and the second column is the parent ID.  
}



\seealso{
 \code{\link{ohg}}, \code{\link{hss}} and \code{\link{rpoh}} 
}
\examples{
set.seed(100)
chr <- list()
sire <- list()
set.seed(1)
chr <- list()
for(i in 1:5)
{
	chr[[i]] <- .simulateHalfsib(numInd = 20, numSNP = 5000, recbound = 1:10)
	sire[[i]] <- ssp(bmh(chr[[i]]), chr[[i]])
	sire[[i]] <- sire[[i]][1,] + sire[[i]][2,]
	sire[[i]][sire[[i]] == 18] <- 9
}

Genotype <- do.call(rbind, chr)
rownames(Genotype) <- 6:(nrow(Genotype) + 5)
sire <- do.call(rbind, sire)
rownames(sire) <- 1:5
Genotype <- rbind(sire, Genotype)
oh <- ohg(Genotype)  # creating the Opposing Homozygote matrix
pogc(oh, 5)
}

\keyword{pedigree}
\keyword{parentage}

