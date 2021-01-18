\name{summaryInc}
\alias{summaryInc}

\title{
Summarize spatial inconsistencies. 
}
\description{
A function to summarize the details of spatial inconsistencies detected with \code{\link{searchInc}}.
}
\usage{
summaryInc(indices, neighborhoods)
}

\arguments{
  \item{indices}{
An integer or integer vector giving the indices of loci (bins) in the corresponding Hi-C matrix for which a summary should be generated 
}
  \item{neighborhoods}{
The list of spatial inconsistencies detected by \code{\link{searchInc}}.\cr
'neighborhoods' is a list where neighborhoods[i] gives the neighborhood of i, if a spatial inconsistency was detected for i, and NULL otherwise. Here, the neighborhood of i is the sub-graph of immediate neighbors in trans, an \code{\link{igraph}} object with attributes to label nodes by their name (index), fullName (chromosome name and genomic coordinates), membership (the connected component to which they belong, where i is always assigned with 0) and chr (chromosome name).
}
}

\value{
A data frame with the following columns, specifying the details for each node in the given detected spatial inconsistencies: chr (chromosome name of node) component (the connected component membership of the node), index (node index), fullName (chromosome name and genomic coordinates of the node), splitChr (chromosome name of the node for which the spatial inconsistency was detected), splitIndex (index of the node for which the spatial inconsistency was detected).
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).}
\author{
Yoli Shavit
}
\seealso{
\code{\link{searchInc}}\cr
\code{\link{plotInc}}\cr
}

\examples{
  data(spatialInc)
  summaryInc(1, spatialInc) # empty since no spatial incosistency was detected
  summaryInc(167, spatialInc)
  summaryInc(167:173, spatialInc)
}
