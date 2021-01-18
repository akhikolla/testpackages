\name{LCT.haplotypes}
\docType{data}
\alias{LCT.haplotypes}
\alias{LCT.hap}
\alias{LCT.sample}
\alias{LCT.snps}

\title{LCT haplotypes data set}

\description{
  These data contain the haplotype matrix \code{LCT.hap} of the 5008 individuals from the 1000 Genomes data for a ~300kb segment containing the Lactase gene.  
  Information about individuals (sex, population and super population) is present in \code{LCT.sample}, and information about snps is available in \code{LCT.snps}.
}

\usage{ data(LCT.haplotypes) }
 
\format{
  Three data objects are present in \code{LCT.haplotypes}:
  \describe{
    \item{\code{LCT.hap}}{ A matrix of haplotypes }
    \item{\code{LCT.sample}}{ A data frame with information on individuals (sex, population, super.population)}
    \item{\code{LCT.snps}}{ A data frame with information on snps (chr, id, dist, pos, A1, A2) }
  }
}

\source{ Data were obtained from the 1000 Genomes Project. }

\seealso{ \code{\link{LCT.matrix}} }

\references{ McVean et al, 2012, \emph{An integrated map of genetic variation from 1,092 human genomes}, Nature \bold{491, 56-65} doi:10.1038/nature11632 }


