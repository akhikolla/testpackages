\name{LCT.matrix}
\docType{data}
\alias{LCT.matrix}
\alias{LCT.matrix.bed}
\alias{LCT.matrix.fam}
\alias{LCT.matrix.pop1000G}

\title{LCT genotypes matrix}

\description{
  These data contain the genotype matrix corresponding to haplotypes present in \code{LCT.haplotypes} from the 1000 Genomes data for a ~300kb segment containing the Lactase gene.  
  Information about individuals is present in \code{LCT.matrix.fam}, and information about population (population and super population) is present in \code{LCT.matrix.pop1000G}, in a format needed to generate a bedmatrix.  
  \code{LCT.snps} from \code{LCT.haplotypes} can be used as the corresponding bim file of this genotypes matrix.  
}

\usage{ data(LCT.matrix) }

\format{
  Three data objects are present in \code{LCT.haplotypes}:
  \describe{
    \item{\code{LCT.matrix.bed}}{ The matrix of genotypes}
    \item{\code{LCT.matrix.fam}}{ The corresponding fam file}
    \item{\code{LCT.matrix.pop1000G}}{ A data frame with population information for individuals (population, superpopulation)}
  }
}

\source{ Data were obtained from the 1000 Genomes Project. }

\seealso{ \code{\link{LCT.haplotypes}} }

\references{ McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632 }

\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)
#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

}
