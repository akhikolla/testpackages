\name{plotInc}
\alias{plotInc}

\title{
Plot a spatial inconsistency 
}
\description{
A function to plot a spatial inconsistency for a given loci (bin).  
}
\usage{
plotInc(i, neighborhoods, label = "fullName", size = 20, interactive = F)
}

\arguments{
  \item{i}{
An integer giving the index of the loci (bin) for which the inconsistency graph should be plotted.
}
  \item{neighborhoods}{
The list of spatial inconsistencies detected by \code{\link{searchInc}}.\cr
'neighborhoods' is a list where neighborhoods[i] gives the neighborhood of i, if a spatial inconsistency was detected for i, and NULL otherwise. Here, the neighborhood of i is the sub-graph of immediate neighbors in trans, an \code{\link{igraph}} object with attributes to label nodes by their name (index), fullName (chromosome name and genomic coordinates), membership (the connected component to which they belong, where i is always assigned with 0) and chr (chromosome name).
}
  \item{label}{
A character vector for labeling the nodes, set to "fullName" by default. Can be one of three options: "name"" (node index), "fullName"" (chromosome name and genomic coordinates) and "chr" (chromosome name). If not one of the above options, then nodes are labelled by their "name" attribute. 
}
  \item{size}{
The size of the nodes in the plot, set to 20 by default. See \code{\link{searchInc}} for more details on plotting options.
}
  \item{interactive}{
A Boolean indicating whether the plot should be interactive or not. See \code{\link{searchInc}} for more details on plotting options.
}
}

\value{
If the chosen neighborhood is not NULL, it is be plotted according to user preferences (nodes are colored by their connected component and labelled based on user choice).
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).}
\author{
Yoli Shavit
}

\seealso{
\code{\link{searchInc}}\cr
\code{\link{summaryInc}}\cr
}

\examples{
  data(spatialInc)
  plotInc(1, spatialInc) # no plot since no spatial incosistency was detected
  plotInc(167, spatialInc) 
}


