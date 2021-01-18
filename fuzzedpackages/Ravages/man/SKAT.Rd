\name{SKAT}
\alias{SKAT}
\title{ SKAT test }
\description{ Peforms SKAT on categorial or binary phenotypes }
\usage{
SKAT(x, NullObject, genomic.region = x@snps$genomic.region, 
     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
     get.moments = "size.based", estimation.pvalue = "kurtosis",
     params.sampling, cores = 10, debug = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF, as it was computed in the original SKAT method}
  \item{maf.threshold}{ The MAF above which variants are removed (default is to keep all variants)}
  \item{get.moments}{ How to estimate the moments to compute the p-values among "size.based", "bootstrap", "permutations", or "theoretical" for categorial phenotypes (2 or more groups of individuals). By default "size.based" that will choose the method depending on sample size (see \code{details}) }
  \item{estimation.pvalue}{ Whether to use the skewness ("skewness") or the kurtosis ("kurtosis") for the chi-square approximation }
  \item{params.sampling}{ A list containing the elements "perm.target", "perm.max", "debug". Only needed if \code{get.moments = "boostrap"} or \code{get.moments = "permutations" } }
  \item{cores}{ How many cores to use for moments computation, set at 10 by default. Only needed if \code{get.moments = "theoretical"} }
  \item{debug}{ Whether to return the mean, standard deviation, skewness and kurtosis of the statistics}
  \item{verbose}{ Whether to display information about the function actions }
}
\details{
  For categorial phenotypes, the p-value is calculated using a chi-square approximation based on the statistics' moments. The user has to choose how to compute these moments (argument \code{get.moments}), and which moments to use for the chi-square approximation (argument \code{estimation.pvalue}).  

  The moments can be computed either using a sampling procedure (\code{"permutations"} if there are no covariates, or \code{"bootstrap"} otherwise), or using theoretical moments computed as in Liu et al. 2008 (\code{"theoretical"}).  

  If \code{get.moments = "size.based"}, the sampling procedure will be used for sample sizes lower than 2000, and the theoretical calculations otherwise.  

  To estimate the p-values, etiher the first three moments are used (\code{estimation.pvalue = "skewness"}), or the moments 1, 2 and 4 are used (\code{estimation.pvalue = "kurtosis"}). 

  If \code{get.moments = "theoretical"} and \code{estimation.pvalue = "skewness"}, it corresponds to \code{method = "liu"} in the SKAT package.
  If \code{get.moments = "theoretical"} and \code{estimation.pvalue = "kurtosis"}, it corresponds to \code{method = "liu.mod"} in the SKAT package.
  
  For small samples, p-values estimation is based on sampling and a sequential procedure: permutated statistics are computed and each one is compared to the observed statistics.
  This method requires \code{perm.target} and \code{perm.max} that should be given as a list to \code{params.bootstrap}.
  If \code{params.bootstrap} is not specified, perm.target will be set at 100, perm.max at 5e4.
  The boostrap progam stops when either \code{perm.target} or \code{perm.max} is reached.  
  P-values are then computed using a mixed procedure:   

  if \code{perm.target} is reached, the p-value is computed as : \code{perm.target} divided by the number of permutations used to reach \code{perm.target};   

  if \code{perm.max} is reached, the SKAT small sample procedure is used, and p-values are approximated using a chi-square distributions based on statistics' moments 1, 2 and 4 computed from the permutated values.     

  If \code{NullObject$pheno.type = "continuous"}, the method from Liu et al. will be used to compute the p-value for the continuous phenotype, but \code{estimation.pvalue} can be set at "skewness" or "kurtosis".  

  If \code{debug=TRUE}, more informations about the estimated statistics moments are given.  
}
\value{
A data frame containing for each genomic region: 
  \item{stat}{ The observed statistics}
  \item{p.value}{ The p-value of the test}

If \code{get.moments = "bootstrap"} or \code{get.moments = "permutations"}, additional fields are present:
  \item{p.perm}{ The p-value computed by permutations: number of times permutated is greater than observed statistics divided by the total number of permutations performed}
  \item{p.chi2}{ The p-value computed by the chi-square approximation using the SKAT small sample procedure}

If \code{debug = TRUE}, the mean, standard deviation, skewness and kurtosis are also returned, as well as for the sampling procedure:
  \item{nb.gep}{ The number of times a permutated statistics is equal or greater than the observed statistics \code{stat}}
  \item{nb.eq}{ The number of times a permutated statistics is equal to the observed statistics \code{stat}}
  \item{nb.perms}{ The total number of simulations performed }
}
\seealso{ \code{\link{NullObject.parameters}}, \code{\link{SKAT.theoretical}}, \code{\link{SKAT.bootstrap}}, \code{\link{SKAT.permutations}}}
\references{ Wu et al. 2011, \emph{Rare-variant association testing for sequencing data with the sequence kernel association test}, American Journal of Human Genetics \bold{82-93} doi:10.1016/j.ajhg.2011.05.029;  

             Lee et al. 2012, \emph{Optimal Unified Approach for Rare-Variant Association Testing with Application to Small-Sample Case-Control Whole-Exome Sequencing Studies}, American Journal of Human Genetics, doi:10.1016/j.ajhg.2012.06.007;  

             Liu et al. 2008, \emph{A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables}, Computational Statistics & Data Analysis, doi:10.1016/j.csda.2008.11.025
}

\examples{
#Example on simulated data from Ravages with 
#One group of 50 controls and 
#two groups of 25 cases, each one with a prevalence of 0.01
#with 50% of causal variants, 5 genomic regions are simulated
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov,
                      n.case.groups = 2, select.gene = "R1",
                      GRR.multiplicative.factor=2)

x.sim <- random.bed.matrix(genes.maf = Kryukov, size = c(50, 25, 25),
                           prev = c(0.001, 0.001), GRR.matrix.del = GRR.del,
                           p.causal = 0.5, p.protect = 0, select.gene="R1",
                           same.variant = FALSE,
                           genetic.model = "multiplicative", replicates = 5)
#Null Model
x.sim.H0 <- NullObject.parameters(x.sim@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")

#Run SKAT (here permutations as n<2000 and no covariates)
#Parameters for the sampling procedure: target = 5, max = 100
#Please increase the number of permutations for a more accurate estimation of the p-values
params.sampling = list(perm.target = 5, perm.max = 100)
SKAT(x.sim, x.sim.H0, params.sampling = params.sampling)

#Run SKAT with a random continuous phenotype
#Null Model
x.sim.H0.c <- NullObject.parameters(rnorm(100), RVAT = "SKAT", pheno.type = "continuous")
SKAT(x.sim, x.sim.H0.c)



\donttest{
#Example on 1000Genome data
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

#run SKAT using the 1000 genome EUR populations as "outcome"
#with very few permutations
#Please increase the permutations for a more accurate estimation of the p-values
#Fit Null model with covariate sex
x1.H0.covar <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorial",
                                     data = covar, formula = ~ sex)

#Run SKAT with the covariates: use boostrap as n<2000
SKAT(x1, x1.H0.covar, params.sampling = params.sampling, get.moments = "bootstrap")

#Run SKAT using theoretical moments (discourage here as n<2000) and 1 core
SKAT(x1, x1.H0.covar, get.moments = "theoretical", cores = 1)
}

}

