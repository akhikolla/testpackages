\name{setGeneTable}
\alias{setGeneTable}
\title{
    setGeneTable
}
\description{
    Takes a \code{data.frame} object with columns 'gene','chr','start', and 'end', and creates an SQLite table of the information. Returns a count of the number of genes in the table.
}
\usage{
    setGeneTable(geneInfo,table='gene',db='snplistdb') 
}
\arguments{
  \item{geneInfo}{
    A \code{data.frame} object of gene location info with columns 'gene','chr','start', and 'end'.
}
  \item{table}{
    Name of the SQLite table to be created. Default is 'gene'.
}
  \item{db}{
    Name of the SQLite database in which to create \code{table}. Default is 'snplistdb'.
}
}

\value{
    Count of genes included in \code{table}.
}

\examples{
    geneInfo <- cbind(c('BRCA1','BRCA2'),c(17,13),c(41196312,32889611),c(41277500,32973805))
    colnames(geneInfo) <- c('gene','chr','start','end')
    setGeneTable(as.data.frame(geneInfo))
}
\keyword{ database }