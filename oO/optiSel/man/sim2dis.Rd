\name{sim2dis}
\alias{sim2dis}

\title{Converts a Similarity Matrix into a Dissimilarity Matrix}

\description{
Converts a similarity matrix (e.g. a kinship matrix) into a dissimilarity matrix.
}
\usage{
sim2dis(f, a=4.0, baseF=0.03, method=1)
}
\arguments{
 \item{f}{Similarity matrix.}
 \item{a}{Exponent}
 \item{baseF}{Old inbreeding not measured by \code{f}}
 \item{method}{Either \code{1} or \code{2}}
}

\value{
Dissimilarity matrix D.
}
\details{
This function converts a similarity matrix \code{f} with values between 0 and 1 (e.g. a kinship matrix) into a dissimilarity matrix.
At first, the similarity is adjusted as

\code{f <- baseF + (1-baseF)*f}.

Then, for Method 1, the dissimilarity between individuals \code{i} and \code{j} is computed as

  \code{Dij = (-log(fij))^a},

whereas for Method 2, the dissimilarity is computed as

\code{Dij = sqrt((fii+fjj)/2-fij)^a}.

Although Method 2 may provide lower stress values in some cases, Method 1 has the advantage that the area reflects the diversity of a population more reasonable.
}
\examples{
data(map)
dir   <- system.file("extdata", package = "optiSel")
files <- file.path(dir, paste("Chr", unique(map$Chr), ".phased", sep=""))
f     <- segIBD(files, map, minSNP=15, minL=1.0)
D     <- sim2dis(f, 4)

## Multidimensional scaling of animals:
\dontrun{
data(Cattle)
library("smacof")
color <- c(Angler="red", Rotbunt="green", Fleckvieh="blue", Holstein="black")
col   <- color[as.character(Cattle$Breed)]
Res   <- smacofSym(D, itmax = 5000, eps = 1e-08)
plot(Res$conf, pch=18, col=col, main="Multidimensional Scaling", cex=0.5)
mtext(paste("segIBD Stress1 = ", round(Res$stress,3)))
}


}
\author{Robin Wellmann}

