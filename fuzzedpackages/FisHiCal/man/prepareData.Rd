\name{prepareData}
\alias{prepareData}
\title{
Match FISH distances and Hi-C frequencies for calibration
}
\description{
A function to match FISH distances and Hi-C frequencies based on matching probe/bin coordinates. }
\usage{
prepareData(fish, fishCoord, hic, hicCoord)
}

\arguments{
  \item{fish}{
A N*N matrix giving the pairwise FISH distances between the probes defined by fishCoord. Off-diagonal zero values stand for missing distances. 
}
  \item{fishCoord}{
The genomic coordinates of the N FISH probes. A data frame with the following columns: chr (chromosome name), start (start position), end (end position). Additional columns are permitted but will not be used.
Note that chromosome names should match between hicCoord and fishCoord.
}
  \item{hic}{
A M*M matrix giving the pairwise Hi-C contact frequencies between the genomic bins defined by hicCoord. Off-diagonal zero values stand for missing frequencies.
}
  \item{hicCoord}{
The genomic coordinates of the M Hi-C bins. A data frame with the following mandatory columns: chr (chromosome name), start (start position), end (end position). Additional columns are permitted but will not be used.
Note that chromosome names should match between hicCoord and fishCoord.
}
}
\details{
This function calls \code{\link{findMatchingIndices}} in order to find matching Hi-C bins for FISH probes.
When several FISH probes are mapped to the same  bin i, the FISH distances for this bin with another bin j are not unique.  In these cases the minimal non-zero FISH distance between i and j is taken as representative, in order to generate a unique match and since Hi-C is likely to be biased towards shorter distances. 
}
\value{
A data frame of matching FISH distances and Hi-C frequencies, sorted in increasing order by distance value. First column is named distances, giving the FISH distances values. Second column is named frequencies, giving the matching Hi-C frequencies values. 
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\author{
Yoli Shavit
}
\seealso{
\code{\link{prepareCalib}}\cr
\code{\link{calibrate}}
}
