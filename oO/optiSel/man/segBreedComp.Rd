\name{segBreedComp}
\alias{segBreedComp}

\title{
Calculates the Segment-Based Breed Composition of Individuals
}

\description{
Calculates the \bold{seg}ment based \bold{Breed} \bold{Comp}osition: 
For every individual the breed composition is estimated, including the genetic contribution from native ancestors.
}
\usage{
segBreedComp(Native, map, unitP="Mb")
}

 

\arguments{
\item{Native}{This parameter is either

(1) Mx(2N) logical matrix, with \code{TRUE}, if the segment containing the SNP is considered native, and \code{FALSE} otherwise. The row names are the marker names, and the non-unique column names are the IDs of the individuals. The matrix is typically computed from component \code{freq} of the output from function \link{haplofreq}.

or

(2) Mx(2N) character matrix, with components being the first characters of the names of the breeds in which the respective segment has maximum frequency. Segments considered native are coded as \code{'1'}. The row names are the marker names, and the non-unique column names are the IDs of the individuals. The matrix is typically component \code{match} from the output of function \link{haplofreq}.

or

(3)  Vector with file names. The files contain for every SNP and for each haplotype 1 if the segment containing the SNP is considered native. Otherwise it is the first letter of the name of the breed in which the segment has maximum frequency. These files are typically created by function \link{haplofreq}.
There is one file per chromosome and file names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.nat"}.
}
\item{map}{
Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in Mega base pairs \code{'Mb'}, and the position in centimorgan \code{'cM'}. The markers must be in the same order as in  in \code{Native}.
}
\item{unitP}{The unit for measuring the proportion of the genome included in native segments.
  Possible units are the number of marker SNPs included in shared segments (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the total length of the shared segments in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
}
\value{
Data frame with the number of rows being the number of individuals. The columns are
\item{Indiv}{IDs of the individuals,}
\item{native}{Genetic contributions from native ancestors,}
\item{...}{Contributions from other breeds.}
}
\details{
For every individual the breed composition is computed, including the genetic contribution from native ancestors (native contribution). The native contribution is the proportion of the genome belonging to segments whose frequency is smaller than a predefined value in all other breeds. 

Additionally, for each introgressed breed,
the proportion of the genome of each individual is computed that is non-native and has maximum frequency in the respective breed (not if option (1) is used).
}
\examples{
data(map)
data(Cattle)
dir     <- system.file("extdata", package = "optiSel")
GTfiles <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
Haplo   <- haplofreq(GTfiles, Cattle, map, thisBreed="Angler", minSNP=20, minL=1.0)
Comp    <- segBreedComp(Haplo$freq<0.01, map)
mean(Comp$native)
#[1] 0.3853432

Comp <- segBreedComp(Haplo$match, map)
apply(Comp[, -1], 2, mean)

## Reading native segments from files:
\dontrun{
wdir <- file.path(tempdir(), "HaplotypeEval")
file <- haplofreq(GTfiles, Cattle, map, thisBreed="Angler", minSNP=20, 
           minL=1.0, ubFreq=0.01, what="match", w.dir=wdir)
Comp <- segBreedComp(file$match, map)
head(Comp)

apply(Comp[, -1], 2, mean)
#    native          F          H          R 
#0.38534317 0.05503451 0.25986508 0.29975724 

#unlink(wdir, recursive = TRUE)
}
}
\author{Robin Wellmann}

