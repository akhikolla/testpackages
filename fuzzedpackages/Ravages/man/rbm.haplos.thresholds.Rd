\name{rbm.haplos.thresholds}
\alias{rbm.haplos.thresholds}
\title{
  Simulation of genetic data based on haplotypes and a libaility model
}
\description{
  Simulates genetic data with respect to allele frequency spectrum and linkage disequilibrium pattern observed on given haplotype data under a libaility model
}
\usage{
  rbm.haplos.thresholds(haplos, weights = -0.4*log10(colMeans(haplos)), 
                        maf.threshold = 0.01, nb.causal, p.protect = 0, 
                        h2, prev, normal.approx = TRUE, size, 
                        replicates, rep.by.causal, verbose = TRUE)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{weights}{ A vector of weights for each variant to compute the burden used in the liability model. By default, wieghts=-0.4*log10(MAF)}
  \item{maf.threshold}{ The maf threshold to consider a rare variant  (set at 0.01 by default), variants with a MAF upper this threshold will have a weight of 0}
  \item{nb.causal}{ The number of causal variants}
  \item{p.protect}{ The proportion of protective variants among causal variants}
  \item{h2}{ The variance explained by the gene}
  \item{prev}{ A vector with the prevalence in each group of individuals}
  \item{normal.approx}{ TRUE/FALSE: whether to use the normal approximation to compute thresholds. Set at TRUE by default}
  \item{size}{ The sizes of each group of individuals}
  \item{replicates}{ The number of simulations to perform}
  \item{rep.by.causal}{ The number of time causal variants will be sampled}
  \item{verbose}{ Whether to display information about the function actions }
}
\value{
  \item{x}{ A bed matrix with simulated genotypes}
}
\details{
  \code{nb.causal}, \code{p.protect}, \code{h2} and \code{prev} should be vectors of length corresponding to the number of groups to simulate. If they are of size 1, values will be duplicated. 

  All monomorphic variants and variants with a MAF higher than \code{maf.threshold} will have a weight of 0. Causal variants are sampled among variants having weights greater than 0.  

  A liability model is built on haplotypes' burden computed on sampled causal variants using each variant's \code{weights}, and adjusted on the desired \code{h2}. Thresholds from this liability are then chosen to respect the given \code{prev} (from a standard normal distribution if \code{normal.approx=TRUE}, or using a distribution from 1e6 sampled burdens if \code{normal.approx=FALSE}). Please be carreful when using the normal approximation with high \code{h2} values or low \code{prev} values. 
  Haplotypes' probabilities in each group of individuals are then computed and two haplotypes are then sampled for each individual based on these probabilities.  

  To simulate a group of controls, \code{prev} needs to be set at 1, regardless of the other arguments.  

  N \code{replicates} will be performed, and to gain in computation time, the same causal variants can be used for multiple replicates as different haplotypes will be sampled for each individual. \code{rep.by.causal} indicates the number of replicates to perform for each set of causal variants.
  To ensure a variability in the simulations, we yet recommend to resample causal variants a few times when many replicates are to be performed. 
  For example, if 1000 replicates are to be performed, we recommend to resample causal variants 20 times.  

  The phenotype will be stored in \code{@ped$pheno}, and the simulation number is \code{@snps$genomic.region}.
}
\examples{
\donttest{
  #Load LCT dataset for haplotype matrix
  data(LCT.haplotypes)
  #LCT gene in the EUR population
  LCT.gene.hap <- LCT.hap[which(LCT.sample$super.population=="EUR"), 
                          which(LCT.snps$pos>=136545410 & LCT.snps$pos<=136594750)]

  #Simulation of 100 controls, and two groups of 50 cases with 30 causal variants
  #and with the second group having half h2 and twice the prevalence
  #compared to the first one
  #5 replicates are performed and causal variants are sampled once
  x <- rbm.haplos.thresholds(haplos=LCT.gene.hap, maf.threshold = 0.01, nb.causal=30,
                             p.protect=0, h2=c(0.01, 0.01, 0.02), prev=c(1, 0.01, 0.005),
                             size=c(100, 50, 50), replicates = 5, rep.by.causal = 5)
}
}
