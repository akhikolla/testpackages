\name{segN}
\alias{segN}

\title{Calculates Probabilities of Alleles to belong to Native Segments}

\description{
\bold{Seg}ment based probability of alleles to be \bold{N}ative:
 For each pair of individuals the probability is computed that two SNPs taken at random position from randomly chosen haplotypes both belong to native segments.
}
\usage{
segN(Native, map, unitP="Mb", keep=NULL, cores=1)
}

\arguments{
\item{Native}{This parameter is either

(1) Mx(2N) indicator matrix, with 1, if the segment containing the SNP is considered native, and 0 otherwise. The row names are the marker names, and the non-unique column names are the IDs of the individuals. The matrix is typically computed from the output of function \link{haplofreq}.

or

(2)  Vector with file names. The files contain for every SNP and for each haplotype from this breed 1 if the segment containing the SNP is considered native. These files are typically created by function \link{haplofreq}.
There is one file per chromosome and file names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.nat"}.
}
\item{map}{
Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in Mega base pairs \code{'Mb'}, and the position in centiMorgan \code{'cM'}. The markers must be in the same order as in \code{Native}.
}
\item{unitP}{The unit for measuring the proportion of the genome included in native segments.
  Possible units are the number of marker SNPs included in shared segments (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the total length of the shared segments in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
 \item{keep}{Vector with IDs of individuals (from this breed) for which the probabilities are to be computed. By default, they will be computed for all individuals included in \code{Native}.}
  \item{cores}{Number of cores to be used for parallel processing of chromosomes. By default one core is used. For \code{cores=NA} the number of cores will be chosen automatically. Using more than one core increases execution time if the function is already fast.}
}
\value{
\code{NxN} matrix with \code{N} being the number of genotyped individuals from this breed (which are also included in vector \code{keep}).
}
\details{
For each pair of individuals the probability is computed that two SNPs taken at random position from randomly chosen haplotypes both belong to native segments. That is, they are not introgressed from other breeds.
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package = "optiSel")
files <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
Freq  <- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="others", minSNP=20)$freq
fN   <- segN(Freq<0.01, map)
mean(fN)
#[1] 0.15418

\dontrun{
fN   <- segN(Freq<0.01, map, cores=NA)
mean(fN)
#[1] 0.15418
}


## using files:
\dontrun{
wdir   <- file.path(tempdir(),"HaplotypeEval")
chr    <- unique(map$Chr)
GTfile <- file.path( dir, paste("Chr", chr, ".phased",     sep=""))
files  <- haplofreq(GTfile, Cattle, map, thisBreed="Angler", w.dir=wdir)

fN     <- segN(files$match, map)
mean(fN)
#[1] 0.15418

fN     <- segN(files$match, map, cores=NA)
mean(fN)
#[1] 0.15418

#unlink(wdir, recursive = TRUE)
}


}
\author{Robin Wellmann}

