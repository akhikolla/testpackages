\name{dils-package}
\alias{dils-package}
\docType{package}
\title{
Data-Informed Link Strength. Combine multiple-relationship networks into a single weighted network.
}
\description{
Combine multiple-relationship networks into a single weighted network.  The approach is similar to factor analysys in the that contribution from each constituent network varies so as to maximize the information gleaned from the multimetwork.  This implementation uses Principal Component Analysis calculated using 'prcomp' with bootstrap subsampling.
}
\details{
\tabular{ll}{
Package: \tab dils\cr
Type: \tab Package\cr
Version: \tab 0.8\cr
Date: \tab 2013-10-27\cr
License: \tab MIT + file LICENSE\cr
}
Start with a table (data.frame, tab-delimited file, database) where each row/record represents a link between two nodes (a dyad) in a directed or undirected network and each column represents a different relationship between the two nodes, ie. each column is a network. DILS combines these columns/networks into a single network that is a weighted sum of the constituent networks.  The resulting DILS network uses information from all of the constituent networks and contains more information than any of the constituent networks.  The output is a data.frame of DILS scores for each dyad, therefore is a single network ready for analysis using \pkg{igraph} or other social network analysis (SNA) tools.

Workflow synthesizing networks might typically look like this:
\enumerate{
  \item Start with several networks in igraph, adjacency list, or edgelist form.
  \item Is necessary, use \code{\link{EdgelistFromIgraph}} or \code{\link{EdgelistFromAdjacency}} to convert igraph and adjacency list networks to edgelist form.
  \item Use \code{\link{MergeEdgelists}} to combine the individual network datasets into a single dataset.
  \item Use \code{\link{GenerateDilsNetwork}} to synthesize the networks in the merged data set into a single weighted network.
  \item Use \code{\link{IgraphFromEdgelist}} or \code{\link{AdjacencyFromEdgelist}} to convert the edgelist output to the desired output.
  \item Use \code{\link{RelativeNetworkInformation}} on input networks and DILS network to see if/how much the information content of the DILS network exceeds the information content of the input networks.
}

Workflow for imputing edges for a binary network might typically look like this:
\enumerate{
  \item Start with a binary network as an adjacency matrix (for an igraph use \code{get.adjacency}).
  \item Use \code{\link{RelationStrengthSimilarity}} to Calculate RSS scores for each dyad.
  \item Use \code{\link{RssSuggestedNetwork}} on the original network and the RelationStrengthSimilarity output to get a new suggested network with more edges.
}
    
}
\author{
Stephen R. Haptonstahl <srh@haptonstahl.org>
}
\references{
"Discovering Missing Links in Networks Using Similarity Measures", 
Hung-Hsuan Chen, Liang Gou, Xiaolong (Luke) Zhang, C. Lee Giles. 2012.
}
\keyword{ network }
