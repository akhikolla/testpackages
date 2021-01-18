\name{CAST}
\alias{CAST}
\title{ Cohort Allelic Sum Test}
\description{ Calculates the CAST genetic score }
\usage{
CAST(x, genomic.region = x@snps$genomic.region, maf.threshold = 0.5)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{maf.threshold}{ The MAF used for the definition of a rare variant, set at 0.5 by default, i.e. all variants are kept}
}
\value{
A matrix containing the CAST genetic score with one row per individual and one column per \code{genomic.region}
}
\seealso{ \code{\link{WSS}}, \code{\link{burden.weighted.matrix}}, \code{\link{burden.mlogit}} }
\references{ Morgenthaler S and Thilly WG. \emph{A strategy to discover genes that carry multi-allelic or mono-allelic risk for common diseases: a cohort allelic sums test (CAST).} Mutat Res. 2007}
\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

# Group variants within known genes
x <- set.genomic.region(x)

# Filter variants with maf (computed on whole sample) < 0.025
# keeping only genomic region with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

# Compute burden score CAST
score.CAST <- CAST(x1, maf.threshold=0.025)
}
