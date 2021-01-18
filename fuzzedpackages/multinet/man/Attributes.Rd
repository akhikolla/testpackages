\name{multinet.attributes}
\alias{multinet.attributes}
\alias{add_attributes_ml}
\alias{attributes_ml}
\alias{set_values_ml}
\alias{get_values_ml}
\title{
Managing attributes
}
\description{
These functions are used to assign and retrieve values to/from actors, vertices.
}
\usage{
add_attributes_ml(n, attributes, type="string", target="actor",
	layer="", layer1="", layer2="")
attributes_ml(n, target="actor")
get_values_ml(n, attribute, actors=character(0),
	vertices =character(0), edges=character(0))
set_values_ml(n, attribute, actors=character(0),
	vertices=character(0), edges=character(0), values)
}
\arguments{
\item{n}{A multilayer network.}
\item{attributes}{Name(s) of the attributes to be created.}
\item{target}{Can be "actor" (attributes attached to actors), "vertex" (attributes attached to vertices) or "edge" (attributes attached to edges). Layer attributes are not available in this version.}
\item{type}{Can be "string" or "numeric".}
\item{layer}{This can be specified only for targets "vertex" (so that the attribute exists only for the vertices in that layer) or "edge" (in which case the attribute applies to intra-layer edges in that layer).}
\item{layer1}{This can be specified only for target "edge", together with layer2, so that the attribute applies to inter-layer edges from layer1 to layer2. If \code{layer1} and \code{layer2} are specified, the parameter \code{layer} should not be used.}
\item{layer2}{See layer1.}
\item{attribute}{The name of the attribute to be updated.}
\item{actors}{A vector of actor names. If this is specified, layers, vertices and edges should not.}
\item{vertices}{A dataframe of vertices to be updated. The first column specifies actor names, the second layer names. If this is specified, actors, layers and edges should not.}
\item{edges}{A dataframe containing the vertices to be connected. The four columns must contain, in this order: actor1 name, layer1 name, actor2 name, layer2 name. If this is specified, actors, layers and vertices should not.}
\item{values}{A vector of values to be set for the corresponding actors, vertices or edges.}
}
\seealso{\link{multinet.properties}, \link{multinet.edge_directionality}}
\value{
\code{attributes_ml} returns a data frame with columns: "name", and "type". If vertex attributes are listed, an additional "layer" column is used. If edge attributes are listed, two columns "layer1" and "layer2" are included.
\code{get_values_ml} returns a data frame with the values for the requested objects.
}
\examples{
net <- ml_aucs()
attributes_ml(net)
# actor attributes, of string type (default)
add_attributes_ml(net,c("name","surname"))
# a numeric attribute associated to the layers (not available in this version)
# add_attributes_ml(net,"num vertices",type="numeric",target="layer")
# attributes for vertices on the facebook layer
add_attributes_ml(net,"username",type="string",target="vertex",layer="facebook")
# attributes for edges on the work layer
add_attributes_ml(net,"strength",type="numeric",target="edge",layer="work")
# listing the attributes
attributes_ml(net)
# attributes_ml(net,"layer") # not available in this version
attributes_ml(net,"vertex")
attributes_ml(net,"edge")
# setting some values for the newly created attributes
set_values_ml(net,"name",actors=c("U54","U139"),values=c("John","Johanna"))
e <- data.frame(
    c("U139","U139"),
    c("work","work"),
    c("U71","U97"),
    c("work","work"))
set_values_ml(net,"strength",edges=e,values=.47)
# getting the values back
get_values_ml(net,"name",actors=c("U139"))
get_values_ml(net,"strength",edges=e)
# setting attributes based on network properties: create a "degree"
# attribute and set its value to the degree of each actor
actors_ml(net) -> a
layers_ml(net) -> l
degree_ml(net,actors=a,layers=l,mode="all") -> d
add_attributes_ml(net,target="actor",type="numeric",attributes="degree")
set_values_ml(net,attribute="degree",actors=a,values=d)
get_values_ml(net,attribute="degree",actors="U54")
# select actors based on attribute values (e.g., with degree greater than 40)
get_values_ml(net,attribute="degree",actors=a) -> degrees
a[degrees>40]
# list all the attributes again
attributes_ml(net)
}
