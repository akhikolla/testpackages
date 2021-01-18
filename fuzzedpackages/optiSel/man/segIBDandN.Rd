\name{segIBDandN}
\alias{segIBDandN}

\title{Calculates Probabilities that Alleles belong to a Shared Native Segment}

\description{
Calculates the \bold{seg}ment based probability of alleles to be \bold{IBD} (identical by descent) \bold{and} \bold{N}ative: 
 For each pair of individuals the probability is computed that two SNPs taken at random position from randomly chosen haplotypes belong to a shared segment and are native.
}
\usage{
segIBDandN(files, Native, map, minSNP=20, unitP="Mb", minL=1.0, 
   unitL="Mb", a=0.0, keep=NULL, skip=NA, cskip=NA, cores=1)
}

 

\arguments{
  \item{files}{Vector with names of the phased marker files, one file for each chromosome. The required format is described under \code{Details}. File names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}.
  }
\item{Native}{This parameter is either

(1) Mx(2N) indicator matrix, with 1, if the segment containing the SNP is considered native, and 0 otherwise. The row names are the marker names, and the non-unique column names are the IDs of the individuals. The matrix is typically computed from the output of function \link{haplofreq}.

or

(2)  Vector with file names. The files contain for every SNP and for each haplotype from this breed 1 if the segment containing the SNP is considered native. These files are typically created by function \link{haplofreq}.
There is one file per chromosome and file names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.nat"}.

}
\item{map}{
Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in mega base pairs \code{'Mb'}, and the position in centimorgan \code{'cM'}. The markers must be in the same order as in \code{files} and in \code{Native}.
}
\item{minSNP}{Minimum number of marker SNPs included in a segment.}
  \item{unitP}{The unit for measuring the proportion of the genome included in shared segments.
  Possible units are the number of marker SNPs included in shared segments (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the total length of the shared segments in centiMorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.
  }
\item{minL}{Minimum length of a segment in \code{unitL} (e.g. in cM or Mb).}
\item{unitL}{The unit for measuring the length of a segment. Possible units are the number of marker SNPs included in the segment (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the genetic  distances between the first and the last marker in centiMorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
  \item{a}{The Function providing the weighting factor for each segment is w(x)=x*x/(a+x*x). The parameter of the function is the length of the segment in \code{unitL}. The default value \code{a=0.0} implies no weighting, whereas \code{a>0.0} implies that old inbreeding has less influence on the result than new inbreeding.}
  \item{keep}{Vector with IDs of individuals (from this breed) for which the probabilities are to be computed. By default, they will be computed for all individuals included in \code{Native}.}
     \item{skip}{Take line \code{skip+1} of the genotype files as the line with column names. By default, the number is determined automatically.}
 \item{cskip}{Take column \code{cskip+1} of the genotype files as the first column with genotypes. By default, the number is determined automatically.}
  \item{cores}{Number of cores to be used for parallel processing of chromosomes. By default one core is used. For \code{cores=NA} the number of cores will be chosen automatically. Using more than one core increases execution time if the function is already fast.}
}
\value{
\code{NxN} matrix with \code{N} being the number of individuals from this breed included in all files (and in parameter \code{keep}).
}
\details{
For each pair of individuals the probability is computed that two SNPs taken at random position from randomly chosen haplotypes belong to a shared segment and are native. That is, they are not introgressed from other breeds.

\bold{Genotype file format}: Each file containing phased genotypes has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. The first \code{cskip} columns are ignored. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces. The name of each file must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}.
}
\examples{
data(map)
data(Cattle)
dir    <- system.file("extdata", package = "optiSel")
GTfile <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
Freq   <- haplofreq(GTfile, Cattle, map, thisBreed="Angler", refBreeds="others", minSNP=20)$freq

fIBDN  <- segIBDandN(GTfile, Freq<0.01, map=map, minSNP=20)
mean(fIBDN)
#[1] 0.01032261

\dontrun{
fIBDN  <- segIBDandN(GTfile, Freq<0.01, map=map, minSNP=20, cores=NA)
mean(fIBDN)
#[1] 0.01032261
}

## using files:
\dontrun{
wdir   <- file.path(tempdir(),"HaplotypeEval")
chr    <- unique(map$Chr)
GTfile <- file.path( dir, paste("Chr", chr, ".phased",     sep=""))
file   <- haplofreq(GTfile, Cattle, map, thisBreed="Angler", minSNP=20, ubFreq=0.01, w.dir=wdir)

fIBDN  <- segIBDandN(GTfile, file$match, map=map, minSNP=20)
mean(fIBDN)
#[1] 0.01032261

fIBDN  <- segIBDandN(GTfile, file$match, map=map, minSNP=20, cores=NA)
mean(fIBDN)
#[1] 0.01032261


#unlink(wdir, recursive = TRUE)
}

}
\author{Robin Wellmann}

