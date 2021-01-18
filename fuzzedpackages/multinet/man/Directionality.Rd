\name{multinet.edge_directionality}
\alias{multinet.edge_directionality}
\alias{set_directed_ml}
\alias{is_directed_ml}
\title{
Controlling edge directionality
}
\description{
Functions to get and set the edge directionality of one or more pairs of layers (that is, the directionality of edges connecting nodes in those layers).
}
\usage{
set_directed_ml(n, directionalities)
is_directed_ml(n, layers1 = character(0), layers2 = character(0))
}
\arguments{
\item{n}{A multilayer network.}
\item{directionalities}{A dataframe with three columns where each row contains a pair of layers (l1,l2) and 0 or 1 (indicating resp. undirected and directed edges). Directionality is automatically set for both (l1,l2) and (l2,l1).}
\item{layers1}{The layer(s) from where the edges start. If \code{layers1} is not provided, all layers are considered.}
\item{layers2}{The layer(s) where the edges end. If an empty list of layers is passed (default), the ending layers are set as equal to those in parameter \code{layers1}.}
}
\value{
\code{is_directed_ml} returns a data frame where each row contains the name of two layers and the corresponding type of edges (directed/undirected).
}
\seealso{\link{multinet.properties}, \link{multinet.attributes}}
\examples{
net <- ml_empty()
# Adding some layers, one directed and one undirected
add_layers_ml(net,c("l1","l2"),c(TRUE,FALSE))
# Setting the directionality of inter-layer edges
layers = c("l1","l2")
dir <- data.frame(layers,layers,c(0,1))
set_directed_ml(net,dir)
# retrieving all directionalities
dir <- is_directed_ml(net)
# copying directionalities to a new network
net2 <- ml_empty()
add_layers_ml(net2,c("l1","l2"))
set_directed_ml(net2,dir)
}
