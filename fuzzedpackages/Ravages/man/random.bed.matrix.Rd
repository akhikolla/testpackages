\name{random.bed.matrix}
\alias{random.bed.matrix}
\title{
  Simulation of genetic data using GRR values
}
\description{
  Generates a simulated bed.matrix with genotypes for cases and controls based on GRR values
}
\usage{
random.bed.matrix(genes.maf = Kryukov, size, prev, replicates, 
                  GRR.matrix.del, GRR.matrix.pro, 
                  p.causal = 0.5, p.protect = 0, same.variant = FALSE, 
                  genetic.model=c("general", "multiplicative", 
                                  "dominant", "recessive"), select.gene)
}
\arguments{
  \item{genes.maf}{ A dataframe containing at least the MAF in the general population (column maf) for variants with their associated gene (column gene), by default the file \code{Kryukov} is used}
  \item{size}{ A vector containing the size of each group (the first one being the control group)}
  \item{prev}{ A vector containing the prevalence of each group of cases}
  \item{replicates}{ The number of simulations to perform }
  \item{GRR.matrix.del}{ A list containing the GRR matrix associated to the heterozygous genotype compared to the homozygous reference genotype as if all variants are deleterious. An additional GRR matrix associated to the homozygous for the alternate allele is needed if \code{genetic.genetic.model="general"}}
  \item{GRR.matrix.pro}{ The same argument as \code{GRR.matrix.del} but for protective variants}
  \item{p.causal}{ The proportion of causal variants in cases}
  \item{p.protect}{ The proportion of protective variants in cases among causal variants}
  \item{same.variant}{ TRUE/FALSE: whether the causal variants are the same in the different groups of cases}
  \item{genetic.model}{ The genetic model of the disease }
  \item{select.gene}{ Which gene to choose from \code{genes.maf$gene} if multiple genes are present. If missing, only the first level is kept.}
}
\value{
  A bed.matrix with as much columns (variants) as \code{replicates}*number of variants.
  The field \code{x@snps$genomic.region} contains the replicate number and the field \code{x@ped$pheno} contrains the group of each individual, "0" being the controls group.
}
\details{
  The genetic model of the disease needs to be specified in this function.

  If \code{genetic.model="general"}, there is no link between the GRR for the heterozygous genotype and the GRR for the homozygous alternative genotype. 
  Therefore, the user has to give two matrices of GRR, one for the heterozygous genotype, the other for the homozygous alternative genotype.

  If \code{genetic.model="multiplicative"}, we assume that the the GRR for the homozygous alternative genotype is the square of the GRR for the heterozygous genotype. 

  If \code{genetic.model="dominant"}, we assume that the GRR for the heterozygous genotype and the GRR for the homozygous alternative genotype are equal.

  If \code{genetic.model="recessive"}, we assume that the GRR for the heterozygous genotype is equal to 1: the GRR given is the one associated to the homozygous alternative genotype.

  \code{GRR.matrix.del} contains GRR values as if all variants are deleterious. These values will be used only for the proportion \code{p.causal} of variants that will be sampled as causal.  

  The files \code{Kryukov} or \code{GnomADgenes} available with the package Ravages can be used as the argument \code{genes.maf}.

  If \code{GRR.matrix.del} (or \code{GRR.matrix.pro}) has been generated using the function \code{GRR.matrix}, the arguments \code{genes.maf} and \code{select.gene} should have
  the same value as in \code{GRR.matrix}.
}
\seealso{ \code{\link{GRR.matrix}}, \code{\link{Kryukov}}, \code{\link{GnomADgenes}} }
\examples{
#GRR values calculated with the SKAT formula
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, 
                      n.case.groups = 2, select.gene = "R1",
                      GRR.multiplicative.factor=2)
                              
#Simulation of one group of 1,000 controls and two groups of 500 cases, 
#each one with a prevalence of 0.001
#with 50% of causal variants, 5 genomic regions are simulated.
x <- random.bed.matrix(genes.maf = Kryukov, size = c(1000, 500, 500), 
                       prev = c(0.001, 0.001), GRR.matrix.del = GRR.del, 
                       p.causal = 0.5, p.protect = 0, select.gene="R1",
                       same.variant = FALSE, 
                       genetic.model = "multiplicative", replicates = 5)
}
