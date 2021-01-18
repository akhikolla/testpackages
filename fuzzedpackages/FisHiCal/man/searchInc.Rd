\name{searchInc}
\alias{searchInc}
\title{
Search for spatial inconsistencies
}
\description{
A function to test for spatial inconsistencies in a calibrated Hi-C distance matrix.
}
\usage{
searchInc(calibHiC, hicCoord)
}
\arguments{
  \item{calibHiC}{
A M*M Hi-C matrix,providing pairwise calibrated distances between M genomic bins. Zero off-diagonal values represent discarded/missing information. This matrix could be prepared with the function \code{\link{calibrate}}.
}
  \item{hicCoord}{
The genomic coordinates of the M Hi-C bins. A data frame with the following mandatory columns: chr (chromosome name), start (start position), end (end position). Additional columns are permitted but will not be used.
}
}
\details{
Denote delta(i,j) the calibrated Hi-C distances, then delta(i,j) define a weighted undirected graph G\{V, E\}, where V is the set of loci (bins) and E is the set of edges: \{(vi,vj)| delta(i,j) > 0, i != j\} with weights delta(i,j). We distinguish here between immediate neighbors from the same chromosome (cis) and from different chromosomes (trans) and detect a spatial inconsistency for a node v in G, if the sub-graph G' of all (immediate) trans neighbors of v is not connected. Further identifying the connected components in G' can highlight the cause of inconsistency and the underlying spatial division.  
}
\value{
A list with the following values:
\item{neighborhoods}{A list of length M, where neighborhoods[i] gives the neighborhood of i, if a spatial inconsistency was detected for i, and NULL otherwise. Here, the neighborhood of i is the sub-graph of immediate neighbors in trans, an \code{\link{igraph}} object with attributes to label nodes by their name (index in hicCoord), fullName (chromosome and genomic coordinates), membership (the connected component to which they belong, where i is always assigned with 0) and chr (chromosome name).}
\item{indices}{The indices of loci (bins) with a spatial inconsistency.}
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).}
\author{
Yoli Shavit
}
\seealso{
\code{\link{prepareCalib}}\cr
\code{\link{calibrate}}\cr
\code{\link{plotInc}}\cr
\code{\link{summaryInc}}\cr
}

