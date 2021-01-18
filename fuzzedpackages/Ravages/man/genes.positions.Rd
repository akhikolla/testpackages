\name{genes.positions}
\docType{data}
\alias{genes.positions}
\alias{genes.b37}
\alias{genes.b38}

\title{Genes positions}

\description{
  Positions of human genes. These data were downloaded from Biomart on the Ensembl website with the GRCh37 and GRCh38 versions. Only genes present in GnomAD were kept.

  Data contain the \code{Chr}, the \code{Start} position, the \code{End} position and the \code{Gene_Name} of all the genes in chromosomes 1 to 22 representing 19375 and 18278 genes in the two GRCh versions respectively.
}
  
\usage{data(genes.b37)
data(genes.b38) }

\format{
  The data contain one dataframe with four columns:
  \describe{
    \item{\code{Chr}}{The chromosome of the gene}
    \item{\code{Start}}{The start position of the gene}
    \item{\code{End}}{The end position of the gene}
    \item{\code{Gene_Name}}{The name of the gene}
  }
}  

\source{The data were obtained from the Ensembl website.}
\references{ RJ Kinsella et al, 2011, \emph{Ensembl BioMarts: a hub for data retrieval across taxonomic space}, Database. doi:10.1093/database/bar030;   

AD Yates et al, 2020, \emph{Ensembl 2020}, Nucleic Acide Research. doi:10.1093/nar/gkz966}
\seealso{ \code{\link{set.genomic.region}} }
