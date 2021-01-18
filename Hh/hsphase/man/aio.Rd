\name{aio}
\alias{aio}
\title{
All-in-one Phasing
}
\description{
Phasing of a half-sib family group.
}
\usage{
aio(genotypeMatrix, bmh_forwardVectorSize = 30, bmh_excludeFP = TRUE,
bmh_nsap = 3, output = "phase")
}
\arguments{
  \item{genotypeMatrix}{
{\code{matrix} half-sib genotypes  (one half-sib per row, with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)
}

}
\item{bmh_forwardVectorSize}{
  \code{integer} number of heterozygous sites used to validate recombination events or check for genotyping errors}


\item{bmh_excludeFP}{
  \code{logical} excludes SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{bmh_nsap}{
  \code{integer} number of SNP per block to validate recombinations}

\item{output}{
  \code{character} if equal to the \code{phase} the 'aio' will only return the phasing results 

}
}
\details{
This function calls the \code{\link{bmh}},  \code{\link{ssp}} and  \code{\link{phf}} functions. 
}
\value{
Returns a list of matrices. The first element (phasedHalfsibs) is a matrix with two rows (phased haplotypes) per individual (first paternal and second maternal). Data in format 0 (A), 1 (B) and 9 (unphased or missing).
The second (sireHaplotype) and third (blockStructure) elements are the same as the output of  \code{\link{ssp}} and  \code{\link{bmh}}.
}

\note{
Only this function needs to be called to phase a half-sib family. The genotype's matrix must contain individuals from only one half-sib family and one ordered chromosome.
}

\seealso{
 \code{\link{bmh}},  \code{\link{ssp}} and  \code{\link{phf}}
}
\examples{
genotype <- matrix(c(       # Define a Half-sib Genotype Matrix
  2,1,0,                    # Individual 1
  2,0,0,                    # Individual 2
  0,0,2                     # Individual 3
  ), byrow = TRUE, ncol = 3)    # There are 3 individulas with three SNPs
  
aio(genotype)               # The genotypes must include only one half-sib family and one chromosome
}

\keyword{snp}
\keyword{phase}
\keyword{haplotype}
\keyword{inference}
\keyword{reconstruction}
