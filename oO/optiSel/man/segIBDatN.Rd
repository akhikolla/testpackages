\name{segIBDatN}
\alias{segIBDatN}

\title{Segment-Based Kinship at Native Alleles.}

\description{
Calculates the kinship at native alleles, which is the segment based probability of native alleles to be IBD.}

\usage{
segIBDatN(files, phen, map, thisBreed, refBreeds="others", ubFreq=0.01, minSNP=20,
  unitP="Mb", minL=1.0, unitL="Mb", a=0.0, keep=NULL, lowMem=TRUE, 
  skip=NA, cskip=NA, cores=1)
}

\arguments{

 \item{files}{This can be a character vector with names of the phased marker files, one file for each chromosome. 
 Alternatively \code{files} can be a list with the following components:

a) \code{hap.thisBreed}: A character vector with names of the phased marker files for the individuals from \code{thisBreed}, one file for each chromosome.
 
b) \code{hap.refBreeds}: A character vector with names of the phased marker files for the individuals from the reference breeds (\code{refBreeds}), one file for each chromosome. If this component is missing, then it is assumed that the haplotypes of these animals are also included in \code{hap.thisBreed}.

c) \code{match}: If present, a character vector with file names containing the origin of the marker alleles. The files are typically created with function \code{haplofreq}. If this vector is missing, then the default method is used to estimate the origins. 

File names must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}. The required format of the marker files is described under \code{Details}.
}

\item{phen}{Data frame containing the ID (column \code{"Indiv"}),  breed name (column \code{"Breed"}), and sex (column \code{Sex}) of each individual.}
\item{map}{Data frame providing the marker map with columns including marker name \code{'Name'}, chromosome number \code{'Chr'}, and possibly the position on the chromosome in Mega base pairs \code{'Mb'}, and the position in centimorgan \code{'cM'}. (The position in base pairs could result in an integer overflow). The order of the markers must bethe same as in the files.}
\item{thisBreed}{Breed name: Results will be computed for individuals from \code{thisBreed}.}
\item{refBreeds}{Vector containing names of genotyped breeds. A segment is considered native if its frequency is smaller than \code{ubFreq}  in all \code{refBreeds}. The default \code{"others"} means that all genotyped breeds except \code{thisBreed} are considered.}
\item{ubFreq}{A segment is considered native if its frequency is smaller than \code{ubFreq} in all reference breeds.}
\item{minSNP}{Minimum number of marker SNPs included in a segment.}
\item{unitP}{The unit for measuring the proportion of the genome included in native segments.
  Possible units are the number of marker SNPs included in shared segments (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the total length of the shared segments in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
\item{minL}{Minimum length of a segment in \code{unitL} (e.g. in cM).}
\item{unitL}{The unit for measuring the length of a segment. Possible units are the number of marker SNPs included in the segment (\code{'SNP'}), the number of Mega base pairs (\code{'Mb'}), and the genetic  distances between the first and the last marker in centimorgan (\code{'cM'}). In the last two cases the map must include columns with the respective names.}
\item{a}{The function providing the weighting factor for each segment is w(x)=x*x/(a+x*x). The parameter of the function is the length of the segment in \code{unitL}. The default value \code{a=0.0} implies no weighting, whereas \code{a>0.0} implies that old inbreeding has less influence on the result than new inbreeding.}
 \item{keep}{
 Subset of the IDs of the individuals from data frame \code{phen} (including individuals from other breeds) or a logical vector indicating the animals in data frame \code{phen} that should be used. By default all individuals included in \code{phen} will be used.}
\item{lowMem}{If \code{lowMem=TRUE} then temporary files will be created and deleted.}
     \item{skip}{Take line \code{skip+1} of the genotype files as the row with column names. By default, the number is determined automatically.}
 \item{cskip}{Take column \code{cskip+1} of the genotype files as the first column with genotypes. By default, the number is determined automatically.}
  \item{cores}{Number of cores to be used for parallel processing of chromosomes. By default one core is used. For \code{cores=NA} the number of cores will be chosen automatically. Using more than one core increases execution time if the function is already fast.}
}
\value{
A list of class \code{ratioFun} including components:
\item{Q1}{matrix with \code{Q1[i,j]}  = Probability that two alleles chosen from individuals i and j are IBD and are native.}
\item{Q2}{matrix with \code{Q2[i,j]} = Probability that two alleles chosen from individuals i and j are both native.}
\item{d1}{The value by which the probability that two alleles chosen from the offspring are IBD and native increases due to genetic drift.}
\item{d2}{The value by which the probability that two alleles chosen from the offspring are native increases due to genetic drift.}
\item{id}{IDs of the individuals for which the probabilites have been computed.}
\item{mean}{Mean kinship at native alleles of the specified individuals. Note that \code{1-mean} is the genetic diversity at native segments of the specified individuals from \code{thisBreed}}

}

\details{
Calculates a list containing matrices needed to compute segment based kinships at native alleles, defined as the conditional probability that two randomly chosen alleles are IBD, given that both originate from native ancestors. An allele is considered to originate from a native ancesor if the segment containing the allele has low frequency in all reference breeds.

The kinship at native alleles between individuals i and j is \code{Q1[i,j]/Q2[i,j]}.

The mean kinship at native alleles in the offspring is \code{(x'Q1x+d1)/(x'Q2x+d2)}, where \code{x} is the vector with genetic contributions of the selection candidates.

\bold{Genotype file format}: Each file containing phased genotypes has a header and no row names. Cells are separated by blank spaces. The number of rows is equal to the number of markers from the respective chromosome and the markers are in the same order as in the \code{map}. The first \code{cskip} columns are ignored. The remaining columns contain genotypes of individuals written as two alleles separated by a character, e.g. A/B, 0/1, A|B, A B, or 0 1. The same two symbols must be used for all markers. Column names are the IDs of the individuals. If the blank space is used as separator then the ID of each individual should repeated in the header to get a regular delimited file. The columns to be skipped and the individual IDs must have no white spaces. The name of each file must contain the chromosome name as specified in the \code{map} in the form \code{"ChrNAME."}, e.g. \code{"Breed2.Chr1.phased"}.
}
\examples{
data(map)
data(Cattle)
dir     <- system.file("extdata", package = "optiSel")
files   <- paste(dir, "/Chr", 1:2, ".phased", sep="")
sKinatN <- segIBDatN(files, Cattle, map, thisBreed="Angler", 
                   ubFreq=0.01, minL=1.0, lowMem=FALSE)
               
## Mean kinship at native segments:
sKinatN$mean
#[1] 0.06695171

## Note that this can not be computed as mean(sKinatN$of).


## Results for individuals:
sKinatN$of <- sKinatN$Q1/sKinatN$Q2
sKinatN$of["Angler1","Angler5"]
#[1] 0.4394066

## Use temporary files to reduce working memory:

sKinatN <- segIBDatN(files, Cattle, map, thisBreed="Angler", ubFreq=0.01,  minL=1.0)
               
## Mean kinship at native segments:
sKinatN$mean
#[1] 0.06695171


}

\author{Robin Wellmann}

