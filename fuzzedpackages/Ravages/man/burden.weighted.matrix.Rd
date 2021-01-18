\name{burden.weighted.matrix}
\alias{burden.weighted.matrix}
\title{ Score matrix for burden tests}
\description{
Computes the score matrix for burden tests based on variants' weights
}
\usage{
burden.weighted.matrix(x, weights, genomic.region = x@snps$genomic.region)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{weights}{ A vector containing the weight of each variant}
  \item{genomic.region}{ A factorcontaining the genomic region of each variant }
}
\value{
A matrix containing the computed genetic score with one row per individual and one column per \code{genomic.region}.
}
\details{
For variant \emph{i} and individual \emph{j}, the genetic score will be computed as weight of variant \emph{i} * number of minor alleles for individual \emph{j}.
}
\seealso{ \code{\link{CAST}}, \code{\link{WSS}}, \code{\link{burden.mlogit}} }
\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

# Group variants within known genes
x <- set.genomic.region(x)

# Filter variants with maf (computed on whole sample) < 0.025
# keeping only genomic region with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

#Compute burden score with weights = 1-maf
score.burden <- burden.weighted.matrix(x1, weights=1-x1@snps$maf)
}

