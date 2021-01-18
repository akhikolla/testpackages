\name{WSS}
\alias{WSS}
\title{ WSS genetic score}
\description{
Caluclates the WSS genetic score
}
\usage{
WSS(x, genomic.region = x@snps$genomic.region)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genomic.region}{ A factor containing the genomic region of each variant }
}
\value{
A matrix containing the WSS genetic score with one row per individual and one column per \code{genomic.region}
}
\seealso{ \code{\link{CAST}}, \code{\link{burden.weighted.matrix}}, \code{\link{burden.mlogit}} }
\references{ Madsen E and Browning S. \emph{A Groupwise Association Test for Rare Mutations Using a Weighted Sum Statistic.} PLoS Genet. 2009 }
\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

# Group variants within known genes
x <- set.genomic.region(x)

# Filter variants with maf (computed on whole sample) < 0.025
# keeping only genomic region with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

# Compute burden score WSS
score.WSS <- WSS(x1)
}
