\name{phf}
\alias{phf}

\title{
Half-Sib Family Phasing
}
\description{
Phases the half-sib family by using the blocking structure and imputed sire matrices.
}
\usage{
phf(GenotypeMatrix, blockMatrix, sirePhasedMatrix)
}

\arguments{
  \item{GenotypeMatrix}{\code{matrix} half-sib genotypes  (one half-sib per row, with SNP ordered by mapping position in the columns. Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)
}
  \item{blockMatrix}{
\code{matrix} blocking structure (output of \code{\link{bmh}})
}
  \item{sirePhasedMatrix}{\code{matrix}  imputed sire (output of \code{\link{ssp}})
}
}
\value{
Returns a matrix that contains the phased parental haplotypes of the half-sibs. It uses 1, 0 and 9 for A, B and missing.
}
\note{
The genotype matrix must only contain individuals from one half-sib family and one ordered chromosome.
This function is used by the \code{\link{aio}} function for complete phasing of a half-sib group.}

\seealso{
\code{\link{aio}}
}
\examples{
genotype <- matrix(c(
  2,1,0,
  2,0,0,
  0,0,2), byrow = TRUE, ncol = 3)
block <- bmh(genotype)
phf(genotype, block, ssp(block, genotype))
}


\keyword{snp}
\keyword{phase}
\keyword{haplotype}
\keyword{inference}
\keyword{reconstruction}
