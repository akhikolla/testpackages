\name{getBioMartData}
\alias{getBioMartData}
\title{
    getBioMartData
}
\description{
    A function leveraging the \code{biomaRt} package to retrieve gene chromosome and start and end positions from Ensembl.
}
\usage{
    getBioMartData(genes,verbose=FALSE,...) 
}
\arguments{
  \item{genes}{
    A vector of gene names matching \code{hgnc_symbol} in the Ensembl database.
}
  \item{verbose}{
    A Boolean indicating whether to output the funcitons progress in terms of the dimensions of the \code{data.frame} being constructed. Default is FALSE.
}
  \item{...}{
    Additional arguments passed on to the internal call to \code{biomaRt::useMart(...)}. If no such arguments are provided, \code{useMart("ensembl", dataset="hsapiens_gene_ensembl")} is run by default.
}
}
\value{
  A \code{data.frame} object with columns 'gene','chr','start', and 'end', suitable for input to the \code{setGeneTable} function.
}
\note{
At the time of package release, the BioMart community portal is temporarily unavailable. See www.biomart.org for updated status or more information. To access alternative hosts, pass additional arguments to the internal call to \code{biomaRt::useMart(...)}, as in the second example below.
}
\references{
  Durinck S., Spellman P.T., Birney E. and Huber W. (2009) Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt, \emph{Nature Protocols}, \bold{4}, 1184--1191.
}
\seealso{
    \code{\link{setGeneTable}}
}
\examples{
\dontrun{
getBioMartData(c("BRCA1","BRCA2"))
getBioMartData(c("BRCA1","BRCA2"), 
               host="www.ensembl.org", 
               biomart="ENSEMBL_MART_ENSEMBL", 
               dataset="hsapiens_gene_ensembl")
}
}
\keyword{ database }

