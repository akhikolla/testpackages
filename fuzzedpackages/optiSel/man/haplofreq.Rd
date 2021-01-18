\name{haplofreq}
\alias{haplofreq}

\title{Evaluates the Occurrence of Haplotype Segments in Particular Breeds}

\description{
For each haplotype from \code{thisBreed} and every SNP the occurence of the haplotype segment containing the SNP in a set of reference breeds is evaluated. The maximum frequency each segment has in one of these reference breeds is computed, and the breed in which the segment has maximum frequency is identified. Results are either returned in a list or saved to files.
}
\usage{
haplofreq(files, phen, map, thisBreed, refBreeds="others", minSNP=20, minL=1.0, 
  unitL="Mb", ubFreq=0.01, keep=NULL, skip=NA, cskip=NA, w.dir=NA, 
  what=c("freq", "match"), cores=1)
}

\arguments{
 \item{files}{Either a character vector with file names, or a list containing character vectors with file names. The files contain phased genotypes,  one file for each chromosome. File names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}. The required format of the marker files is described under \code{Details}.
 
 If \code{file} is a character vector then, genotypes of all animals must be in the same files. Alternatively, \code{files} can be a list with the following two components:
%\tabular{ll}{

 \code{hap.thisBreed}: Character vector with names of the phased marker files for the individuals from \code{thisBreed}, one file for each chromosome.
 
  \code{hap.refBreeds}: Character vector with names of the phased marker files for the individuals from the reference breeds (\code{refBreeds}), one file for each chromosome. If this component is missing, then it is assumed that the haplotypes of these animals are also included in \code{hap.thisBreed}.
  
%}
}
\item{phen}{Data frame containing the ID (column \code{"Indiv"}) and the breed name (column \code{"Breed"}) of each genotyped individual.}
\item{map}{Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in mega base pairs \code{'Mb'}, and the position in centimorgan \code{'cM'}. The order of the markers must be the same as in the files \code{files}. Marker names must have no white spaces.}
\item{thisBreed}{Name of a breed from column \code{Breed} in \code{phen}: The occurence of each haplotype segment from this breed in the reference breeds will be evaluated.}
\item{refBreeds}{Vector with names of breeds from column \code{Breed} in \code{phen}. These breeds are used as reference breeds.  The occurence of haplotype segments in these breeds will be evaluated. By default, all breeds in \code{phen}, except \code{thisBreed} are used as reference breeds. In contrast, for \code{refBreeds="all"}, all genotyped breeds are used as reference breeds.}
\item{minSNP}{Minimum number of marker SNPs included in a segment.}
\item{minL}{Minimum length of a segment in \code{unitL} (e.g. in cM or Mb).}
\item{unitL}{The unit for measuring the length of a segment. Possible units are the number of marker SNPs included in the segment (\code{'SNP'}), the number of mega base pairs (\code{'Mb'}), and the genetic  distances between the first and the last marker in centiMorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
\item{ubFreq}{If a haplotype segment has frequency smaller than \code{ubFreq} in all reference breeds then the breed name is replaced by \code{'1'}, which indicates that the segment is native.}
 \item{keep}{Subset of the IDs of the individuals from data frame \code{phen}, or a logical vector indicating the animals in data frame \code{phen} that should be used.
 The default \code{keep=NULL} means that all individuals included in \code{phen} will be considered.}
 \item{skip}{Take line \code{skip+1} of the files as the line with column names. By default, the number is determined automatically.}
 \item{cskip}{Take column \code{cskip+1} of the files as the first column with genotypes. By default, the number is determined automatically.}
 \item{w.dir}{Output file directory. Writing results to files has the advantage that much less working memory is required. By default, no files are created. The function returns only the file names if files are created. }
 \item{what}{For \code{what="freq"}, the maximum frequency each haplotype segment has in the reference breeds will be computed. For \code{what="match"}, the name of the reference breed in which the segment has maximum frequency will be determined. By default, the frequencies and the breed names both are determined.}
 \item{cores}{Number of cores to be used for parallel processing of chromosomes. By default one core is used. For \code{cores=NA} the number of cores will be chosen automatically. Using more than one core increases execution time if the function is already fast.}
 }
\value{
If \code{w.dir=NA} then a list is returned. The list may have the following components:
\item{freq}{Mx(2N) - matrix containing for every SNP and for each of the 2N haplotypes from \code{thisBreed} the maximum frequency the segment containing the SNP has in a the reference breeds.}
\item{match}{Mx(2N) - matrix containing for every SNP and for each of the 2N haplotypes from \code{thisBreed} the first letter of the name of the reference breed in which the segment containing the SNP has maximum frequency. Segments with frequencies smaller than \code{ubFreq} in all reference breeds are marked as \code{'1'}, which indicates that the segment is native for \code{thisBreed}.
}
The list has attributes \code{thisBreed}, and \code{map}. 

If \code{w.dir} is the name of a directory, then results are written to files, whereby each file corresponds to one  chromosome, and a data frame with file names is returned.
}
\details{
For each haplotype from \code{thisBreed} and every SNP the occurence of the haplotype segment containing the SNP in a set of reference breeds is evaluated. The maximum frequency each segment has in one of these reference breeds is computed, and the breed in which the segment has maximum frequency is identified. Results are either returned in a list or saved to files.

\bold{Marker file format}: Each marker file containing phased genotypes has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. The first \code{cskip} columns are ignored. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces. 
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package="optiSel")
files <- file.path(dir, paste("Chr", 1:2, ".phased", sep=""))

Freq <- freqlist(
 haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Rotbunt",   minL=2.0),
 haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Holstein",  minL=2.0),
 haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Fleckvieh", minL=2.0)
  )

plot(Freq, ID=1, hap=2, refBreed="Rotbunt")
plot(Freq, ID=1, hap=2, refBreed="Holstein", Chr=1)

\dontrun{
## Test for using multiple cores:

Freq1 <- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Rotbunt", 
                   minL=2.0, cores=NA)$freq
range(Freq[[1]]-Freq1)
#[1] 0 0
}

## Creating output files with allele frequencies and allele origins:
\dontrun{
rdir  <- system.file("extdata", package = "optiSel")
wdir  <- file.path(tempdir(), "HaplotypeEval")
chr   <- unique(map$Chr)
files <- file.path(rdir, paste("Chr", chr, ".phased", sep=""))
wfile <- haplofreq(files, Cattle, map, thisBreed="Angler", minL=2.0, w.dir=wdir)

View(read.table(wfile$match[1],skip=1))
#unlink(wdir, recursive = TRUE)
}

}
\author{Robin Wellmann}

