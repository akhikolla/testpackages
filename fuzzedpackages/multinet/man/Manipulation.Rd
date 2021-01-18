\name{multinet.update}
\alias{multinet.update}
\alias{add_layers_ml}
\alias{add_vertices_ml}
\alias{add_edges_ml}

\alias{add_igraph_layer_ml}

\alias{delete_layers_ml}
\alias{delete_actors_ml}
\alias{delete_vertices_ml}
\alias{delete_edges_ml}

\alias{add_nodes_ml}
\alias{delete_nodes_ml}
\title{
Manipulation of multilayer networks
}
\description{
Functions to add or remove objects in a multilayer network.

The functions \code{add_vertices_ml} and \code{delete_vertices_ml} add/remove the input actors to/from the input layers. Since version 3.1, the actors in the network correspond to the union of all the actors in the various layers (that is, the vertices).

A layer can also be added from an igraph object, where the vertex attribute \code{name} represents the actor name, using the \code{add_igraph_layer_ml} function.

The functions \code{add_nodes_ml} and \code{delete_nodes_ml} are deprecated in the current version of the library. The names vertex/vertices are now
preferentially used over node/nodes.
}
\usage{
add_layers_ml(n, layers, directed=FALSE)
add_vertices_ml(n, vertices)
add_edges_ml(n, edges)

add_igraph_layer_ml(n, g, name)

delete_layers_ml(n, layers)
delete_actors_ml(n, actors)
delete_vertices_ml(n, vertices)
delete_edges_ml(n, edges)
}
\arguments{
\item{n}{A multilayer network.}
\item{layers}{An array of names of layers.}
\item{actors}{An array of names of actors.}
\item{g}{An igraph object with simple edges and a vertex attribute called name storing the actor name corresponding to the vertex.}
\item{name}{Name of the new layer.}
\item{directed}{Determines if the layer(s) is (are) directed or undirected. If multiple layers are specified, directed should be either a single value or an array with as many values as the number of layers.}
\item{vertices}{A dataframe of vertices to be updated or deleted. The first column specifies actor names, the second layer names.}
\item{edges}{A dataframe containing the edges to be connected or deleted. The four columns must contain, in this order: actor1 name, layer1 name, actor2 name, layer2 name.}
}
\seealso{
\link{multinet.properties}, \link{multinet.edge_directionality}
}
\examples{
net <- ml_empty()
# Adding some layers
add_layers_ml(net,"l1")
add_layers_ml(net,c("l2","l3"),c(TRUE,FALSE))
layers_ml(net)
# Adding some vertices (actor A3 is not present in layer l3: no corresponding vertex there)
vertices <- data.frame(
    c("A1","A2","A3","A1","A2","A3"),
    c("l1","l1","l1","l2","l2","l2"))
add_vertices_ml(net,vertices)
vertices <- data.frame(
    c("A1","A2"),
    c("l3","l3"))
add_vertices_ml(net,vertices)
vertices_ml(net)
# Verifying that the actors have been added correctly
num_actors_ml(net)
actors_ml(net)
# We create a data frame specifying two edges:
# A2,l2 -- A3,l1
# A2,l2 -- A3,l2
edges <- data.frame(
    c("A2","A2"),
    c("l2","l2"),
    c("A3","A3"),
    c("l1","l2"))
add_edges_ml(net,edges)
edges_ml(net)

# The following deletes layer 1, and also deletes
# all vertices from "l1" and the edge with an end-point in "l1"
delete_layers_ml(net,"l1")
# The following also deletes the vertices associated to
# "A1" in layers "l2" and "l3"
delete_actors_ml(net,"A1")
# deleting vertex A2,l3 and edge A2,l2 -- A3,l2
delete_vertices_ml(net,data.frame("A2","l3"))
edges <- data.frame("A2","l2","A3","l2")
delete_edges_ml(net,edges)
net
}
