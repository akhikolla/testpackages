\name{read.indiv}
\alias{read.indiv}

\title{Reads Individual IDs from a Genotype File}

\description{
Reads individual IDs from a genotype file.
}
\usage{
read.indiv(file, skip=NA, cskip=NA)
}
\arguments{
 \item{file}{Name of the genotype file.}
  \item{skip}{Take line \code{skip+1} of the genotype files as the row with column names. By default, the number is determined automatically.}
 \item{cskip}{Take column \code{cskip+1} of the genotype files as the first column with genotypes. By default, the number is determined automatically.}
}
\value{
Vector with the IDs of the individuals.
}
\details{
Reading individual IDs from phased marker files. 

\bold{Marker file format}: Each marker file containing phased genotypes has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. The first \code{cskip} columns are ignored. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces. The name of each file must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}.
}
\examples{
data(Cattle)

dir   <- system.file("extdata", package = "optiSel")
file  <- file.path(dir, "Chr1.phased")
ID    <- read.indiv(file)

identical(Cattle$Indiv, ID)
#[1] TRUE


}
\author{Robin Wellmann}

