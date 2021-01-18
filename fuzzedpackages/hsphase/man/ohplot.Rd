\name{ohplot}
\alias{ohplot}
\title{
Opposing Homozygotes Plot
}
\description{
Plot the sorted vectorized matrix of Opposing Homozygotes.
}
\usage{
ohplot(oh, genotype, pedigree, check = FALSE)
}
\arguments{
  \item{oh}{
  \code{integer} Opposing homozygotes matrix (Output of \code{\link{ohg}}) }
  
  \item{genotype}{
\code{matrix} genotype of one chromosome (data should be numeric. Use 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data)
}

 \item{pedigree}{
 \code{matrix}  the pedigree matrix should contain at least two columns, the first column with the half-sib IDs
and the second column with the sires IDs. This argument is optional.
}
  \item{check}{
\code{logical} check the genotype file for the possible errors
}
}
\details{
The \code{cut off} line shows the edge of most different groups. 
}



\seealso{
\code{\link{ohg}} and \code{\link{rpoh}}
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
pedigree <- as.matrix(data.frame(c(1:5, 6:(nrow(Genotype))), 
rep = c(rep(0,5), rep(1:5, rep(20,5)))))
ohplot(oh, Genotype, pedigree, check = TRUE)
}

\keyword{plot}
\keyword{Opposing_Homozygotes}

