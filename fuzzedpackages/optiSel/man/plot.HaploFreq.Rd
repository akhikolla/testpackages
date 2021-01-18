\name{plot.HaploFreq}
\alias{plot.HaploFreq}

\title{Plots Frequencies of Haplotype Segments in Specified Breeds}

\description{
For a particular haplotype from \code{thisBreed} and each marker \code{m} the frequency of the segment containing marker \code{m} in a specified reference breed is plotted.
}
\usage{
\method{plot}{HaploFreq}(x, ID=1, hap=1, refBreed=NULL, Chr=NULL, show.maxFreq=FALSE, ...)
}

\arguments{
\item{x}{This is either an R-Object obtained with function \link{haplofreq} or a list obtained with function \link{freqlist}. 
}
\item{ID}{Either the ID of the animal from this breed to be plotted, or the position of the animal in R-Object \code{x}.}
\item{hap}{Number of the haplotype to be plotted (\code{1} or \code{2})}
\item{refBreed}{Breed name. The frequencies  the haplotype segments have in this reference breed will be plotted. Parameter \code{refBreeds="others"} means that the maximum frequency will be plotted the segments have in other breeds.
}
\item{Chr}{Vector with chromosomes to be plotted. The default means that all chromosomes will be plotted.}
\item{show.maxFreq}{If \code{show.maxFreq=TRUE} then a peak of the grey curve means that a haplotype segment exist in the breed which has high frequency in one of the reference breeds. This frequency is shown. The default is \code{FALSE}.
 }
\item{...}{Arguments to be passed to methods, such as graphical parameters.}
}
\details{
 For a particular haplotype from \code{thisBreed} and each marker \code{m} from chromosomes \code{Chr} the frequency of the segment containing marker \code{m} in reference breed \code{refBreed} is plotted (red line), as well as the maximum frequency the segment has in one of the evaluated breeds (black line), and the maximum frequency a segment from \code{thisBreed} has in one of the evaluated breeds (grey area, if \code{show.maxFreq=TRUE}).
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package="optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")

Freq <- freqlist(
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Rotbunt",   minSNP=20),
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Holstein",  minSNP=20),
  haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Fleckvieh", minSNP=20)
  )

names(Freq)

plot(Freq, ID=1, hap=2, refBreed="Rotbunt")

Freq <- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="others",   minSNP=20)

plot(Freq, ID=1, hap=2)
plot(Freq, ID=1, hap=2, show.maxFreq=TRUE)

Freq <- haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Angler",   minSNP=20)
plot(Freq, ID=1, hap=2)

}
\author{Robin Wellmann}

