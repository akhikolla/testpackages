\name{burden.mlogit}
\alias{burden.mlogit}
\title{
  Logistic or multinomial regression on a genetic score
}
\description{
  Performs a logistical or a non-ordinal multinomial regression on a genetic score
}
\usage{
burden.mlogit(x, NullObject, genomic.region = x@snps$genomic.region, burden,
              maf.threshold = 0.5, get.OR.value = FALSE, alpha = 0.05, cores = 10)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor containg the genomic region of each SNP, \code{x@snps$genomic.region} by default, only needed if \code{burden}="CAST" or \code{burden}="WSS"}
  \item{burden}{"CAST" or "WSS" to directly compute the CAST or the WSS genetic score; or a matrix with one row per individual and one column per \code{genomic.region} if another genetic score is wanted.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score. Set at 0.5 by default}
  \item{get.OR.value}{ TRUE/FALSE: whether to return OR values}
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{cores}{ How many cores to use for moments computation, set at 10 by default }
}
\details{
  This function will return results from the regression of the phenotype on the genetic score for each genomic region.

  If only two groups of individuals are present, a classical logistic regression is performed. 
  If more than two groups of individuals are present, a non-ordinal multinomial regression is performed, 
  comparing each group of individuals to the reference group indicated by the argument \code{ref.level} in \code{NullObject.parameters}. 
  The choice of the reference group won't affect the p-values, but only the Odds Ratios.
  In both types of regression, the p-value is estimated using the Likelihood Ratio test.

  If another genetic score than CAST or WSS is wanted, a matrix with one row per individual and one column per \code{genomic.region} 
  containing this score should be given to \code{burden}. In this situation, no bed matrix \code{x} is needed.
}
\value{
  A dataframe with one row per genomic region and at least two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{get.OR.value=TRUE}, additional columns are present:
  \item{OR}{ The OR value(s) associated to the regression. If there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR}
  \item{l.upper}{ The upper bound of the confidence interval of each OR}
}
\seealso{ \code{\link{NullObject.parameters}}, \code{\link{CAST}}, \code{\link{WSS}}, \code{\link{burden.weighted.matrix}} }
\references{ Bocher O, et al. DOI: 10.1002/gepi.22210. \emph{Rare variant associationtesting for multicategory phenotype.} Genet.Epidemiol. 2019;43:646–656.}
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

#Simulation of a covariate + Sex as a covariate
sex <- x1@ped$sex
set.seed(1) ; u <- runif(nrow(x1))
covar <- cbind(sex, u)
              
#run null model, using the 1000Genome population as "outcome"
#Null model with the covariate sex
x1.H0.covar <- NullObject.parameters(pheno = x1@ped$pop, ref.level = "CEU", 
                                     RVAT = "burden", pheno.type = "categorial",
                                     data = covar, formula = ~ sex)
#WSS test
burden.mlogit(x1, NullObject = x1.H0.covar, burden = "WSS", get.OR.value=TRUE, cores = 1)

}   
