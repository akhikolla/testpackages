\name{Chr1.phased}
\docType{data}
\alias{Chr1.phased}
\title{Phased Cattle Genotypes from Chromosome 1}
\description{
 Phased genotypes of cattle from chromosome 1 (only the first part of the chromosome). Further information on these animals is included in data frame \link{Cattle}.
}
\format{
All functions reading phased genotype data assume that the files are in the following format:

Each file has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. There can be some extra columns on the left hand side containing no genotype data. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should be repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces.

Use function \link{read.indiv} to extract the IDs of the individuals from a genotype file. 
}
\examples{
GTfile <- system.file("extdata/Chr1.phased", package="optiSel")
file.show(GTfile)
GT     <- read.table(GTfile, header=TRUE, skip=2, check.names=FALSE)
GT[1:10,1:5]
}

\keyword{datasets}