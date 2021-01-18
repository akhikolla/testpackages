\name{impute}
\alias{impute}
\title{
Impute of Low Density SNP Marker to High Density (Paternal Strand)}
\description{
Impute the paternal strand from low density to high density utilising high density sire haplotype.}
\usage{
impute(halfsib_genotype_ld, sire_hd, bmh_forwardVectorSize = 30, 
bmh_excludeFP = TRUE, bmh_nsap = 3)
}
\arguments{
\item{halfsib_genotype_ld}{\code{matrix}
half-sib genotypes with low density marker (one half-sib per row, with SNP ordered by mapping position in the columns.
 Data should be numeric. Use 0, 1 and 2 respectively for AA, AB and BB. Use 9 for missing data)}

\item{sire_hd}{\code{matrix}
haplotype of sire (this parameter can be sequence data or any phased sire - the matrix should have rownames which are the sample
 IDs and colnames which are the SNP names)}

\item{bmh_forwardVectorSize}{
  \code{integer} number of heterozygous sites used to validate recombination events or check for genotyping errors}


\item{bmh_excludeFP}{
  \code{logical} exclude SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{bmh_nsap}{
  \code{integer} number of SNPs per block}
}

\value{
Return an imputed half-sib matrix. 
}

\seealso{
 \code{\link{bmh}},  \code{\link{ssp}} and  \code{\link{phf}}
}

\keyword{impute}
\keyword{Low_Density}
\keyword{High_Density}
