\name{burden.continuous}
\alias{burden.continuous}
\title{
  Linear regression on a genetic score
}
\description{
  Performs a linear regression on a genetic score
}
\usage{
  burden.continuous(x, NullObject, genomic.region = x@snps$genomic.region, 
                    burden, maf.threshold = 0.5)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor containg the genomic region of each SNP, \code{x@snps$genomic.region} by default, only needed if \code{burden}="CAST" or \code{burden}="WSS"}
  \item{burden}{"CAST" or "WSS" to directly compute the CAST or the WSS genetic score, or a matrix with one row per individual and one column per \code{genomic.region} if another genetic score is wanted.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score. Set at 0.5 by default}
}
\details{
  This function will return results from the regression of the continuous phenotype on the genetic score for each genomic region.

  If another genetic score than CAST or WSS is wanted, a matrix with one row per individual and one column per \code{genomic.region}
  containing this score should be given to \code{burden}. In this situation, no bed matrix \code{x} is needed.
}
\value{
  A dataframe with one row per genomic region and at least two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
}
\seealso{ \code{\link{CAST}}, \code{\link{WSS}}, \code{\link{burden.weighted.matrix}} }
\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

#Group variants within known genes
x <- set.genomic.region(x)

#Filter of rare variants: only non-monomorphic variants with
#a MAF lower than 2.5%
#keeping only genomic regions with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

#run burden test WSS, using a random continuous variable as phenotype
x1@ped$pheno <- rnorm(nrow(x1))
#Null model
x1.H0 <- NullObject.parameters(pheno = x1@ped$pheno, 
                               RVAT = "burden", pheno.type = "continuous")
burden.continuous(x1, NullObject = x1.H0, burden = "WSS")

}

