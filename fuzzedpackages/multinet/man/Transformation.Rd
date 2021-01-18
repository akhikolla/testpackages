\name{multinet.transformation}
\alias{multinet.transformation}
\alias{flatten_ml}
\alias{project_ml}
\title{
Functions to transform existing layers into new ones.
}
\description{
These functions merge multiple layers into one. The new layer is added to the network. If the input layers are no longer necessary, they must be explicitely erased.

\code{flatten_ml} adds a new layer with the actors in the input layers and an edge between A and B if they are connected in any of the merged layers.

\code{project_ml} adds a new layer with the actors in the first input layer and an edge between A and B if they are connected to the same actor in the second layer.
}
\usage{
flatten_ml(n, new.layer = "flattening", layers = character(0),
  method = "weighted", force.directed = FALSE, all.actors = FALSE)
project_ml(n, new.layer = "projection", layer1, layer2,
method = "clique")
}
\arguments{
\item{n}{A multilayer network.}
\item{new.layer}{Name of the new layer.}
\item{layers}{An array of layers belonging to the network.}
\item{layer1}{Name of a layer belonging to the network.}
\item{layer2}{Name of a layer belonging to the network.}
\item{method}{This argument can take values "weighted" or "or" for \code{flatten_ml} and "clique"
    for \code{project_ml}. "weighted" adds an attribute to the new edges with the number of
    layers where the two actors are connected.}
\item{force.directed}{The new layer is set as directed. If this is false, the new layer is set as directed if at least one of the merged layers is directed.}
\item{all.actors}{If TRUE, then all the actors are included in the new layer, even if they are not present in any of the merged layers.}

}
\seealso{\link{multinet.conversion}}
\references{
Dickison, Magnani, and Rossi, 2016. Multilayer Social Networks. Cambridge University Press. ISBN: 978-1107438750
}
\examples{
net <- ml_aucs()
# A new layer is added to the network, with a flattening of all the other layers
flatten_ml(net, layers = layers_ml(net))
# Bipartite network
from_actor=c("A","B")
to_actor=c("1","1")
from_layer=c("l1","l1")
to_layer=c("l2","l2")
edges = data.frame(from_actor, from_layer, to_actor, to_layer)
n = ml_empty()
add_edges_ml(n, edges)
project_ml(n, layer1 = "l1", layer2="l2")
}
