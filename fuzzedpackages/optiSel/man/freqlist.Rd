\name{freqlist}
\alias{freqlist}

\title{Combines Objects Computed with Function haplofreq() into a List}

\description{
The function combines objects computed with function \link{haplofreq} into a list with class \code{HaploFreq} and adds some attributes.
}
\usage{
freqlist(...)
}

\arguments{
 \item{...}{R-objects computed with function \link{haplofreq}.}
}
\value{
A list with class \code{HaploFreq}
}
\details{
The function combines objects computed with function \link{haplofreq} into a list with class \code{HaploFreq}.
}
\examples{
data(map)
data(Cattle)
dir   <- system.file("extdata", package="optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")

Freq <- freqlist(
 haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Rotbunt",   minL=2.0),
 haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Holstein",  minL=2.0),
 haplofreq(files, Cattle, map, thisBreed="Angler", refBreeds="Fleckvieh", minL=2.0)
  )

#The component names are the reference breeds by default:
names(Freq)

plot(Freq, ID=1, hap=2, refBreed="Rotbunt")

plot(Freq, ID=1, hap=2, refBreed="Holstein", Chr=1)

}
\author{Robin Wellmann}

