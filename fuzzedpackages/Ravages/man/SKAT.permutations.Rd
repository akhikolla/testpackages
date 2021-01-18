\name{SKAT.permutations}
\alias{SKAT.permutations}
\title{ Multi group SKAT test using bootstrap sampling }
\description{ Peforms SKAT on two or more groups of individuals using bootstrap sampling}
\usage{
SKAT.permutations(x, NullObject, genomic.region = x@snps$genomic.region,
                  weights = (1-x@snps$maf)**24, maf.threshold = 0.5,
                  perm.target = 100, perm.max = 5e4, debug = FALSE,
                  estimation.pvalue = "kurtosis")
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF, as it was computed in the original SKAT method}
  \item{maf.threshold}{ The MAF above which variants are removed (default is to keep all variants)}
  \item{perm.target}{ The number of times to exceed the observed statistics. If not reached, \code{perm.max} permutations will be used }
  \item{perm.max}{ The maximum number of permutations to perform to estimate the p-value, will be used if \code{perm.target} is not reached}
  \item{debug}{ Whether to print details about the permutations (mean, standard deviation, skewness, kurtosis), FALSE by default}
  \item{estimation.pvalue}{ Whether to use the skewness ("skewness") or the kurtosis ("kurtosis") for the chi-square approximation }
}
\details{
  P-values estimation is based on permutations sampling and a sequential procedure: permutated statistics are computed and each one is compared to the observed statistics.
  The boostrap progam stops when either \code{perm.target} or \code{perm.max} is reached.
  P-values are then computed using a mixed procedure:

  if \code{perm.target} is reached, the p-value is computed as : \code{perm.target} divided by the number of permutations used to reach \code{perm.target};

  if \code{perm.max} is reached, p-values are approximated using a chi-square distributions based on the first three moments if \code{estimation.pvalue = "skewness"}, or on statistics' moments 1, 2 and 4 if \code{estimation.pvalue = "kurtosis"}.  

  If \code{debug=TRUE}, more informations about the estimated statistics moments are given.  
  
  This function is used by \code{SKAT} when the sample size is smaller than 2000 and no covariates are present.
}
\value{
A data frame containing for each genomic:
  \item{stat}{ The observed statistics}
  \item{p.value}{ \code{p.perm} if \code{perm.target} is reached, \code{p.chi2} if \code{perm.max} is reached. }
  \item{p.perm}{ The p-value computed by permutations: number of times permutated is greater than observed statistics divided by the total number of permutations performed}
  \item{p.chi2}{ The p-value computed by the chi-square approximation using the SKAT small sample procedure}

If \code{debug=TRUE}, other informations are given about the moments estimation:
  \item{nb.gep}{ The number of times a permutated statistics is equal or greater than the observed statistics \code{stat}}
  \item{nb.eq}{ The number of times a permutated statistics is equal to the observed statistics \code{stat}}
  \item{nb.perms}{ The total number of simulations performed }
  \item{mean}{ The mean of the permutated statistics}
  \item{sigma}{ The standard deviation of the permutated statistics}
  \item{skewness}{ The skweness of the permutated statistics}
  \item{kurtosis}{ The kurtosis of the permutated statistics}
}
\seealso{ \code{\link{NullObject.parameters}}, \code{\link{SKAT}} }
\references{ Wu et al. 2011, \emph{Rare-variant association testing for sequencing data with the sequence kernel association test}, American Journal of Human Genetics \bold{82-93} doi:10.1016/j.ajhg.2011.05.029;  

             Lee et al. 2012, \emph{Optimal Unified Approach for Rare-Variant Association Testing with Application to Small-Sample Case-Control Whole-Exome Sequencing Studies}, American Journal of Human Genetics, doi:10.1016/j.ajhg.2012.06.007;
}
\examples{
\donttest{
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
#a MAF lower than 1%
#keeping only genomic regions with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.01, min.nb.snps = 10)

#run SKAT using the 1000 genome EUR populations as "outcome"
#The maximum number of permutations used is 100,
#and the target number is 10, please increase
#both values for a more accurate estimation of the p-values
#Fit Null model
x1.H0 <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorial")
SKAT.permutations(x1, x1.H0, perm.target = 10, perm.max=100)
}
}
