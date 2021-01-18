\name{genotypic.freq}
\alias{genotypic.freq}
\title{
  Genotypic frequencies calculation for data simulations
}
\description{
  Calculates the three genotypic frequencies in the controls group and each group of cases based on MAF in the general population and GRR values
}
\usage{
genotypic.freq(genes.maf = Kryukov, GRR.het, GRR.homo.alt, prev, 
               genetic.model = c("general", "multiplicative", 
                                 "dominant", "recessive"), select.gene)
}
\arguments{
  \item{genes.maf}{ A file containing the MAF in the general population (column maf) for variants with their associated gene (column gene), by default the file \code{Kryukov} is used}
  \item{GRR.het}{ A matrix giving the GRR of the heterozygous genotype compared to the homozygous reference genotype with one row per cases group and one column per variant }
  \item{GRR.homo.alt}{ A matrix giving the GRR of the homozygous alternative genotype compared to the homozygous reference genotype with one row per cases group and one column per variant, only need if \code{genetic.model="general"}}
  \item{prev}{ A vector containing the prevalence of each group of cases}
  \item{genetic.model}{ The genetic model of the disease }
  \item{select.gene}{ Which gene to choose from \code{genes.maf$gene} if multiple genes are present. If missing, only the first level is kept.}
}
\value{
  A matrix of MAF values with one column per variant and one row per group (the first one being the controls group)
}
\details{
  This function is used to simulate genetic data.

  The genetic model of the disease needs to be specified to \code{genetic.model}:

  If \code{genetic.model="general"}, there is no link between the GRR associated 
  to the heterozygous genotype and the GRR associated to the homozygous alternative 
  genotype. Therefore, the user has to give two matrices of GRR, one for each of these genotypes.

  If \code{genetic.model="multiplicative"}, we assume that the GRR associated to the homozygous 
  alternative genotype is the square of the GRR associated to the heterozygous genotype.

  If \code{genetic.model="dominant"}, we assume that the GRR associated to the heterozygous genotype 
  and the GRR associated to the homozygous alternative genotype are equal.

  If \code{genetic.model="recessive"}, we assume that the GRR associated to the heterozygous genotype 
  is equal to 1: the GRR given is the one associated to the homozygous alternative genotype.  

  \code{prev} corresponds to the proportion of each sub-group of cases in the population. 
  It is used only to calculate the MAF in the controls group.  

  The dataframes \code{Kryukov} or \code{GnomADgenes} available with the package Ravages can be used for the argument \code{genes.maf}.
}
\seealso{ \code{\link{GRR.matrix}}, \code{\link{random.bed.matrix}}, \code{\link{GnomADgenes}}, \code{\link{Kryukov}} }
\examples{
#Construction of the GRR matrix using the formula from SKAT
#to compute the GRR (higher weights to rarer variants)
#GRR in the second group are twice as high as in the first group
GRR.del <- GRR.matrix(GRR = "SKAT", GRR.multiplicative.factor=2,
                      select.gene="R1")

#Calculation of frequency in the three groups of individuals 
#under a multilpicative model of the disease
geno.freq.groups <- genotypic.freq(genes.maf = Kryukov, GRR.het = GRR.del, 
                                   prev = c(0.001, 0.001), select.gene="R1", 
                                   genetic.model = "multiplicative")
}
