\name{multinet.actor_measures}
\alias{multinet.actor_measures}
\alias{degree_ml}
\alias{degree_deviation_ml}
\alias{neighborhood_ml}
\alias{xneighborhood_ml}
\alias{connective_redundancy_ml}
\alias{relevance_ml}
\alias{xrelevance_ml}
\title{
Network analysis measures
}
\description{
These functions compute network analysis measures providing a basic description of the actors in the network.
}
\usage{
degree_ml(n, actors = character(0), layers = character(0), mode = "all")
degree_deviation_ml(n, actors = character(0),
  layers = character(0), mode = "all")
neighborhood_ml(n, actors = character(0),layers = character(0), mode = "all")
xneighborhood_ml(n, actors = character(0),layers = character(0), mode = "all")
connective_redundancy_ml(n, actors = character(0),
  layers = character(0), mode = "all")
relevance_ml(n, actors = character(0),layers = character(0), mode = "all")
xrelevance_ml(n, actors = character(0),layers = character(0), mode = "all")
}
\arguments{
\item{n}{A multilayer network.}
\item{actors}{An array of names of actors.}
\item{layers}{An array of names of layers.}
\item{mode}{This argument can take values "in", "out" or "all" to count respectively incoming edges, outgoing edges or both.}
}
\value{
\code{degree_ml} returns the number of edges adjacent to the input actor restricted to the specified layers.
\code{degree_deviation_ml} returns the standard deviation of the degree of an actor on the input layers. An actor with the same
degree on all layers will have deviation 0, while an actor with a lot of neighbors on one layer and only a few on another will have
a high degree deviation, showing an uneven usage of the layers (or layers with different densities).

\code{neighborhood_ml} returns the number of actors adjacent to the input actor restricted to the specified layers. \code{xneighborhood_ml} returns the number of actors adjacent to the input actor restricted to the specified layers and not present in the other layers.

\code{connective_redundancy_ml} returns 1 minus neighborhood divided by degree_

\code{relevance_ml} returns the percentage of neighbors present on the specified layers. \code{xrelevance_ml} returns the percentage of neighbors present on the specified layers and not on others.
}
\references{
\itemize{
\item Berlingerio, Michele, Michele Coscia, Fosca Giannotti, Anna Monreale, and Dino Pedreschi. 2011. "Foundations of Multidimensional Network Analysis." In International Conference on Social Network Analysis and Mining (ASONAM), 485-89. IEEE Computer Society.
\item Magnani, Matteo, and Luca Rossi. 2011. "The ML-Model for Multi-Layer Social Networks." In International conference on Social Network Analysis and Mining (ASONAM), 5-12. IEEE Computer Society.
}
}
\examples{
net <- ml_aucs()
# degrees of all actors, considering edges on all layers
degree_ml(net)
# degree of actors U54 and U3, only considering layers work and coauthor
degree_ml(net,c("U54","U3"),c("work","coauthor"),"in")
# an indication of whether U54 and U3 are selectively active only on some layers
degree_deviation_ml(net,c("U54","U3"))
# co-workers of U3
neighborhood_ml(net,"U3","work")
# co-workers of U3 who are not connected to U3 on other layers
xneighborhood_ml(net,"U3","work")
# percentage of neighbors of U3 who are also co-workers
relevance_ml(net,"U3","work")
# redundancy between work and lunch
connective_redundancy_ml(net,"U3",c("work","lunch"))
# percentage of neighbors of U3 who would no longer
# be neighbors by removing this layer
xrelevance_ml(net,"U3","work")
}
