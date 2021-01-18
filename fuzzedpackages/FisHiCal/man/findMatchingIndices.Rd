\name{findMatchingIndices}
\alias{findMatchingIndices}
\title{
Find matching Hi-C bins for FISH probes, using mid-points
}
\description{
A function for finding matching indices of Hi-C bins for FISH probes, based on the mid-points of their coordinates. Given the start and end coordinates of a set of FISH probes, the function finds the matching Hi-C bins and returns their indices in the given hicCoord data frame. }
\usage{
findMatchingIndices(fishCoord, hicCoord)
}

\arguments{
  \item{fishCoord}{
The genomic coordinates of the FISH probes. A data frame with the following columns: chr (chromosome name), start (start position), end (end position). Additional columns are permitted but will not be used.
Note that chromosome names should match between hicCoord and fishCoord.
}
  \item{hicCoord}{
The genomic coordinates of the Hi-C bins. A data frame with the following mandatory columns: chr (chromosome name), start (start position), end (end position). Additional columns are permitted but will not be used.
Note that chromosome names should match between hicCoord and fishCoord.
}
}
\value{
The indices of matching Hi-C bins in the given hicCoord data frame.  
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\author{
Yoli Shavit
}
\seealso{
\code{\link{prepareData}}\cr
}
