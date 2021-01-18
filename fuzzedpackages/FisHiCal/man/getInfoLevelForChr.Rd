\name{getInfoLevelForChr}
\alias{getInfoLevelForChr}
\title{
Compute the information level for a given chromosome after calibration
}
\description{
A function to compute the information proportions in a a given calibrated Hi-C matrix for a given chromosome. This will be the proportions of non-zero cells, representing information that was not discarded in the calibration step. 
}
\usage{
getInfoLevelForChr(calibHiC, hicCoord, chr)
}
\arguments{
  \item{calibHiC}{
A calibrated Hi-C matrix. This matrix could be prepared with \code{\link{calibrate}}.
}
 \item{hicCoord}{
The genomic coordinates of the Hi-C bins. A data frame with the following mandatory columns: chr (chromosome name), start (start position), end (end position). Additional columns are permitted but will not be used.
}
\item{chr}{
The name of the chromosome of interest (as it appears in hicCoord).
}

}

\details{
The level of information can affect the accuracy of the 3D prediction, which may be lower for sparse calibration matrices (information level < 0.4). }                                                                   

\value{
The information level for the given chromosome in the calibrated Hi-C matrix.
}

\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\author{
Yoli Shavit
}
\seealso{
\code{\link{calibrate}}\cr
\code{\link{lsmacof}}\cr
}



