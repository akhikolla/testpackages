\name{ssp}
\alias{ssp}
\title{
Sire Imputation and Phasing
}
\description{
Infer (impute) and phase sire's genotype based on the block structure matrix (recombination blocks) and homozygous sites of the half-sib genotype matrix.}
\usage{
ssp(blockMatrix, genotypeMatrix)
}

\arguments{
  \item{blockMatrix}{
\code{matrix} block structure (Output of \code{\link{bmh}})
}
  \item{genotypeMatrix}{
\code{matrix} half-sibs genotype
(each row includes the SNP of individuals, 0, 1 and 2 for respectively AA, AB and BB. Use 9 for missing data) }

}

\value{
Returns a matrix (Imputed Sire) with two rows one for each haplotype of the sire (columns are SNP in the order of the genotype matrix). Alleles are coded as 0 (A) and 1 (B). Alleles that could not be imputed are coded as 9.
}



\seealso{
 \code{\link{phf}}, \code{\link{aio}} and \code{\link{imageplot}}
}
\examples{
genotype <- matrix(c(
0,2,1,1,1,
2,0,1,2,2,
2,2,1,0,2,
2,2,1,1,1,
0,0,2,1,0), ncol = 5, byrow = TRUE) # each row contains the SNP of individuals

(result <- ssp(bmh(genotype), genotype)) 
}
\keyword{snp}
\keyword{halfsib}
\keyword{phase}
