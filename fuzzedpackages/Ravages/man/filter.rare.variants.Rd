\name{filter.rare.variants}
\alias{filter.rare.variants}
\title{
  Rare variants filtering
}
\description{
  Filter rare variants based on a MAF threshold and a given number of SNP per genomic region
}
\usage{
filter.rare.variants(x, ref.level, filter=c("whole", "controls", "any"), 
                     maf.threshold=0.01, min.nb.snps, group, genomic.region = NULL)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{ref.level}{ The level corresponding to the controls group, only needed if \code{filter=="controls"} }
  \item{filter}{ On which group the filter will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant, set at 0.01 by default}
  \item{min.nb.snps}{ The minimum number of snps needed to keep a genomic region}
  \item{group}{ A factor indicating the group of each individual, only needed if \code{filter = "controls"} or \code{filter = "any"}. If missing, \code{x@ped$pheno} is taken}
  \item{genomic.region}{ An optional factor containing the genomic region of each variant, only needed if \code{min.nb.snps} is specified and if \code{x@snps$genomic.region} doesn't exist}
}
\value{
  A bed.matrix with filtered variants
}
\details{
  To use this function, a factor 'genomic.region' should be present in the slot \code{x@snps}.  

  If \code{filter="whole"}, only the variants having a MAF lower than the threshold in the entire sample are kept.

  If \code{filter="controls"}, only the variants having a MAF lower than the threshold in the controls group are kept.

  If \code{filter="any"}, only the variants having a MAF lower than the threshold in any of the groups are kept.
}
\examples{
#Import 1000Genome data from region around LCT gene
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Group variants within known genes
x <- set.genomic.region(x)
table(x@snps$genomic.region, useNA="ifany")
                             
#Filter of rare variants: only non-monomorphic variants with 
#a MAF lower than 2.5%
#keeping only genomic regions with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)
table(x1@snps$genomic.region, useNA="ifany")
}
