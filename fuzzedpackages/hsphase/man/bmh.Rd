\name{bmh}
\alias{bmh}
\title{
Block Partitioning
}
\description{
Identifies the block structure (chromosome segments) in the half-sib family that each individual inherited from its sire.
}
\usage{
bmh(GenotypeMatrix, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3)
}

\arguments{
  \item{GenotypeMatrix}{
          \code{matrix} half-sib genotypes  (one half-sib per row, with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)}
  
  \item{forwardVectorSize}{
  \code{integer} number of heterozygous sites used to validate recombination events or check for genotyping errors (50k -> 30, 700k -> 120)}


\item{excludeFP}{
  \code{logical} excludes SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{nsap}{
  \code{integer} number of SNP per block to validate recombinations (50k -> 3, 700k -> 10)}
}
\value{
Returns a matrix of the blocking structure that contains 1s, 2s and 0s. 1s and 2s are the two sire strands. The choice of strand is arbitrary for each chromosome and not consistent across chromosomes. 0s indicate regions of unknown origin.}
\note{
The genotype's matrix must contain individuals from only one half-sib family and one ordered chromosome.}

\seealso{
 \code{\link{ssp}}, \code{\link{phf}}, \code{\link{aio}} and \code{\link{imageplot}}
}
\examples{
genotype <- matrix(c(
0,2,1,1,1,
2,0,1,2,2,
2,2,1,0,2,
2,2,1,1,1,
0,0,2,1,0), ncol = 5, byrow = TRUE) 

(result <- bmh(genotype))

}

\keyword{snp}
\keyword{halfsib}
\keyword{phase}
