\name{multinet-package}
\alias{multinet-package}
\alias{multinet}
\docType{package}
\title{
Multilayer social network analysis and mining
}
\description{
This package defines a class to store multilayer networks and functions to pre-process, analyze and mine them.

With \emph{multilayer social network} we indicate a network where vertices (V) are organized into multiple layers (L) and each node corresponds to an actor (A), where the same actor can be mapped to nodes in different layers. Formally, a multilayer social network as implemented in this package is a graph G = (V, E) where V is a subset of A x L.

In this manual, \emph{\link{multinet.IO}} describes functions to read and write multilayer networks from/to file and the file format. To quickly test some features of the library, some existing multilayer networks are also included (\emph{\link{multinet.predefined}}). A synthetic multilayer network can be generated using the growing models described in \emph{\link{multinet.generation}}.

Updating and getting information about the basic components of a multilayer network (layers, actors, vertices and edges) can be done using the methods described in \emph{\link{multinet.properties}}, \emph{\link{multinet.update}} and \emph{\link{multinet.edge_directionality}}. \emph{\link{multinet.navigation}} shows how to retrieve the neighbors of a node. Attribute values can also be attached to the basic objects in a multilayer network (actors, layers, vertices and edges). Attribute management is described in \emph{\link{multinet.attributes}}.

Each individual layer as well as combination of layers obtained using the data pre-processing (flattening) functions described in \emph{\link{multinet.transformation}} can be analyzed as a single-layer network using the iGraph package, by converting them as shown in \emph{\link{multinet.conversion}}. We can also visualize small networks using the method described in \emph{\link{multinet.plotting}} and the layouts in \emph{\link{multinet.layout}}.

Multilayer network analysis measures are described in \emph{\link{multinet.actor_measures}} (for single-actor, degree-based measures), \emph{\link{multinet.distance}} (for measures based on geodesic distances) and \emph{\link{multinet.layer_comparison}} (to compare different layers).

Communities can be extracted using various clustering algorithms, described in \emph{\link{multinet.communities}}.

Most of the methods provided by this package are described in the book "Multilayer Social Networks". These methods have been proposed by many different authors: extensive references are available in the book, and in the documentation of each function we indicate the main reference we have followed for the implementation. For a few methods developed after the book was published we give specific references to the corresponding literature.
}
\author{
Matteo Magnani \email{matteo.magnani@it.uu.se}
}
\references{
Dickison, Magnani, and Rossi, 2016. Multilayer Social Networks. Cambridge University Press. ISBN: 978-1107438750
}
\keyword{ multilayer network social analysis mining }
