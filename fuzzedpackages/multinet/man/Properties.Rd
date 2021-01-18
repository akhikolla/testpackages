\name{multinet.properties}
\alias{multinet.properties}
\alias{actors_ml}
\alias{layers_ml}
\alias{vertices_ml}
\alias{edges_ml}
\alias{edges_idx_ml}
\alias{num_actors_ml}
\alias{num_layers_ml}
\alias{num_vertices_ml}
\alias{num_edges_ml}

\alias{nodes_ml}
\alias{num_nodes_ml}
\title{
Listing network properties
}
\description{
These functions are used to list basic information about the components of a multilayer network (actors, layers, vertices and edges). 

The functions \code{nodes_ml} and \code{num_nodes_ml} are deprecated in the current version of the library. The names vertex/vertices are now preferentially used over node/nodes.
}
\usage{
layers_ml(n)
actors_ml(n, layers = character(0))
vertices_ml(n, layers = character(0))
edges_ml(n, layers1 = character(0), layers2 = character(0))
edges_idx_ml(n)

num_layers_ml(n)
num_actors_ml(n, layers = character(0))
num_vertices_ml(n, layers = character(0))
num_edges_ml(n, layers1 = character(0), layers2 = character(0))
}
\arguments{
\item{n}{A multilayer network.}
\item{layers}{An array of names of layers belonging to the network. Only the actors/vertices in these layers are returned. If the array is empty, all the vertices in the network are returned. Notice that this may not correspond to the list of actors: there can be actors that are not present in any layer. These would be returned only using the \code{actors_ml} function.}
\item{layers1}{The layer(s) from where the edges to be extracted start. If an empty list of layers is passed (default), all the layers are considered.}
\item{layers2}{The layer(s) where the edges to be extracted end. If an empty list of layers is passed (default), the ending layers are set as equal to those in parameter layer1.}
}
\value{
\code{actors_ml} and \code{layers_ml} return an array of respectively actor and layer names. \code{vertices_ml} returns a data frame where each row contains the name of the actor corresponding to that vertex and the layer of the vertex. \code{edges_ml} returns a data frame where each row contains two actor names (i.e., an edge), the name of the two layers connected by the edge (which can be the same layer if it is an intra-layer edge) and the type of edge (directed/undirected).

    \code{edges_idx_ml} returns the index of the vertex as returned by the \code{vertices_ml} function instead of its name - this is used internally by the plotting function.

The functions num_* compute the number of objects of the requested type.}
\seealso{
\link{multinet.attributes},
\link{multinet.update}, \link{multinet.edge_directionality}
}
\examples{
net <- ml_aucs()
actors_ml(net)
layers_ml(net)
vertices_ml(net)
# only vertices in the "facebook" layer
vertices_ml(net,"facebook")
# all edges
edges_ml(net)
# Only edges inside the "lunch" layer
edges_ml(net,"lunch","lunch")
# Does the same as in the previous line
edges_ml(net,"lunch")
# Returns an empty  data frame, because there are no edges from the
# "lunch" layer to the "facebook" layer
edges_ml(net,"lunch","facebook")

num_actors_ml(net)
num_layers_ml(net)
num_vertices_ml(net)
# Only vertices in the "facebook" layer are counted
num_vertices_ml(net,"facebook")
num_edges_ml(net)
# Only edges inside the "lunch" layer are counted
num_edges_ml(net,"lunch","lunch")
# Does the same as in the previous line
num_edges_ml(net,"lunch")
# Returns 0, because there are no edges from the "lunch" layer to
# the "facebook" layer
num_edges_ml(net,"lunch","facebook")
}
