\name{rbm.haplos.freqs}
\alias{rbm.haplos.freqs}
\title{
  Simulation of genetic data based on haplotypic frequencies
}
\description{
  Simulates genetic data with respect to allele frequency spectrum and linkage disequilibrium pattern observed on given haplotypes and their frequencies
}
\usage{
  rbm.haplos.freqs(haplos, freqs, size, replicates)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{freqs}{ A matrix of haplotypes frequencies in each group of individuals}
  \item{size}{ The sizes of each group of individuals}
  \item{replicates}{ The number of simulations to perform}
}
\value{
  \item{x}{ A bed matrix with simulated genotypes}
}
\details{
  Simulations are performed to respect linkage disequilibrium pattern and allelic frequency spectrum in each group of individuals
  The phenotypic values will be the colnames of \code{freqs} and stored in \code{@ped$pheno}. The simulation number will be in \code{@snps$genomic.region}.
}
\examples{
  #Simulations of 5 groups of individuals with haplotypes frequencies
  #from the 5 EUR populations

  #Load LCT dataset for haplotype matrix
  data(LCT.haplotypes)
  #Haplotypes for the variants in the LCT gene in the EUR population
  LCT.gene.hap <- LCT.hap[which(LCT.sample$super.population=="EUR"),
                         which(LCT.snps$pos>=136545410 & LCT.snps$pos<=136594750)]

  #Individuals from EUR
  LCT.sample.EUR <- subset(LCT.sample, super.population=="EUR")
  #Matrix of haplotypic frequencies
  LCT.freqs <- sapply(unique(LCT.sample.EUR$population), function(z) 
                      ifelse(LCT.sample.EUR$population==z, 
                             1/table(LCT.sample.EUR$population)[z], 0))

  #Simulation of genetic data for five groups of 50 individuals
  x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size=rep(50,5), replicates=5)
}
