\name{segInbreeding}
\alias{segInbreeding}

\title{Calculates Segment Based Inbreeding}

\description{
\bold{Seg}ment based \bold{Inbreeding}: For each individual the probability is computed that the paternal allele and the maternal allele, sampled from random position, belong to a shared segment (i.e. a run of homozygosity, ROH). The arguments are the same as for function \link{segIBD}.
}
\usage{
segInbreeding(files, map, minSNP=20, minL=1.0, unitP="Mb", unitL="Mb", 
   a=0.0, keep=NULL, skip=NA, cskip=NA)
}
\arguments{
  \item{files}{This parameter is either

(1) A vector with names of phased marker files, one file for each chromosome,

or

(2) A list with two components. Each component is a vector with names of phased marker files, one file for each chromosome. Each components corresponds to a different set of individuals. 

File names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}. The required format of the marker files is described under \code{Details}.
  }
  \item{map}{Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in Mega base pairs \code{'Mb'}, and the position in centimorgan \code{'cM'}. (The position in base pairs could result in an integer overflow.) The order of the markers must be the same as in the files.}
  \item{minSNP}{Minimum number of marker SNPs included in a segment.}
  \item{minL}{Minimum length of a segment in \code{unitL} (e.g. in cM or Mb).}
  \item{unitP}{The unit for measuring the proportion of the genome included in shared segments.
  Possible units are the number of marker SNPs included in shared segments (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the total length of the shared segments in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.
  }
  \item{unitL}{The unit for measuring the length of a segment. Possible units are the number of marker SNPs included in the segment (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the genetic  distances between the first and the last marker in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
  \item{a}{The Function providing the weighting factor for each segment is w(x)=x*x/(a+x*x). The parameter of the function is the length of the segment in \code{unitL}. The default value \code{a=0.0} implies no weighting, whereas \code{a>0.0} implies that old inbreeding has less influence on the result than new inbreeding.}
  \item{keep}{If \code{keep} is a vector containing IDs of individuals then inbreeding will be computed only for these individuals. The default \code{keep=NULL} means that inbreeding will be computed for all individuals included in the files.}
   \item{skip}{Take line \code{skip+1} of the files as the row with column names. By default, the number is determined automatically.}
 \item{cskip}{Take column \code{cskip+1} of the files as the first column with genotypes. By default, the number is determined automatically.}
}
\value{
A data frame with column \code{Indiv} containing the individual IDs and column \code{Inbr} containing the inbreeding coefficients.
}
\details{
  For each pair of individuals the probability is computed that two SNPs taken at random position from randomly chosen haplotypes belong to a shared segment.
  
\bold{Genotype file format}: Each file containing phased genotypes has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. The first \code{cskip} columns are ignored. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces. 
}
\references{
  de Cara MAR, Villanueva B, Toro MA, Fernandez J (2013). Using genomic tools to maintain diversity and fitness in conservation programmes. Molecular Ecology. {22}: 6091-6099
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package = "optiSel")
files <- file.path(dir, paste("Chr", 1:2, ".phased", sep=""))
f     <- segInbreeding(files, map, minSNP=20, minL=2.0)

Cattle2 <- merge(Cattle, f, by="Indiv")
tapply(Cattle2$Inbr, Cattle2$Breed, mean)
#    Angler  Fleckvieh   Holstein    Rotbunt 
#0.03842552 0.05169508 0.12431393 0.08386849 

boxplot(Inbr~Breed, data=Cattle2, ylim=c(0,1), main="Segment Based Inbreeding")

fIBD  <- segIBD(files, map, minSNP=20, minL=2.0)
identical(rownames(fIBD), f$Indiv)
#[1] TRUE

range(2*diag(fIBD)-1-f$Inbr)
#[1] -2.220446e-16  2.220446e-16
}
\author{Robin Wellmann}

