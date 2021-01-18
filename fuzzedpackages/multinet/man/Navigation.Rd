\name{multinet.navigation}
\alias{multinet.navigation}
\alias{neighbors_ml}
\alias{xneighbors_ml}
\title{
Functions to extract neighbors of vertices, to navigate the network
}
\description{
These functions return actors who are connected to the input actor through an edge. They can be used to navigate the graph, following paths inside it.
}
\usage{
neighbors_ml(n, actor, layers = character(0), mode = "all")
xneighbors_ml(n, actor, layers = character(0), mode = "all")
}
\arguments{
\item{n}{A multilayer network.}
\item{actor}{An actor name present in the network, whose neighbors are extracted.}
\item{layers}{An array of layers belonging to the network. Only the nodes in these layers are returned. If the array is empty, all the nodes in the network are returned.}
\item{mode}{This argument can take values "in", "out" or "all" to indicate respectively neighbors reachable via incoming edges, via outgoing edges or both.}
}
\value{
\code{neighbors_ml} returns the actors who are connected to the input actor on at least one of the specified layers. \code{xneighbors_ml} (eXclusive neighbors) returns the actors who are connected to the input actor on at least one of the specified layers, and on none of the other layers. Exclusive neighbors are those neighbors that would be lost by removing the input layers.
}
\seealso{
\link{multinet.properties}
}
\references{
Berlingerio, Michele, Michele Coscia, Fosca Giannotti, Anna Monreale, and Dino Pedreschi. 2011. "Foundations of Multidimensional Network Analysis." In International Conference on Social Network Analysis and Mining (ASONAM), 485-89. IEEE Computer Society.
}
\examples{
net <- ml_aucs()
# out-neighbors of U54, that is, all A such that there is an edge ("U54",A)
neigh <- neighbors_ml(net, "U54", mode="out")
# all in-neighbors of U54 on the "work" layer who are not in-neighbors
# in any other layer
xneigh <- xneighbors_ml(net, "U54", "work", mode="in")
# all neighbors (in- and out-) of U54 on the "work" and "lunch" layers
# who are not neighbors in any other layer
xneigh <- xneighbors_ml(net, "U54", c("work","lunch"))
}
