\name{multinet.IO}
\alias{multinet.IO}
\alias{read_ml}
\alias{write_ml}
\title{
Reading and writing multilayer networks from/to file
}
\description{
These functions are used to store a multilayer network to a file or load it from a file.

There are two network formats accepted: multiplex (default) or multilayer. A full multiplex network input file has the following format:

\verb{
 -- comment lines start with two dashes (--)
 #VERSION
 3.0

 #TYPE
 multiplex

 #ACTOR ATTRIBUTES
 AttributeName1,STRING
 AttributeName2,NUMERIC
 -- etc.

 #NODE ATTRIBUTES
 LayerName1,AttributeName1,STRING
 LayerName1,AttributeName2,NUMERIC
 LayerName2,AttributeName3,STRING
 -- etc.

 #EDGE ATTRIBUTES
 -- edge attributes can be defined for specific layers (called local attributes):
 LayerName1,AttributeName,STRING
 LayerName1,AttributeName,NUMERIC
-- or for all layers (called global attributes):
AnotherAttributeName,NUMERIC
 -- etc.

 #LAYERS
 LayerName1,UNDIRECTED
 LayerName2,DIRECTED
 LayerName3,UNDIRECTED,LOOPS
 -- etc.
 -- LOOPS indicates that edges from one vertex to itself (called loops) are allowed on that layer

 #ACTORS
 ActorName1,AttributeValueList...
 ActorName2,AttributeValueList...
 -- etc.

 #VERTICES
 ActorName1,LayerName1,AttributeValueList...
 ActorName1,LayerName2,AttributeValueList...
 -- etc.

 #EDGES
 ActorName1,ActorName2,LayerName1,LocalAttributeValueList,GlobalAttributeValueList...
 -- etc.
 -- the attribute values must be specified in the same order in which they are defined above

 ----------------------------------
}

 If the #LAYERS section is empty, all edges are created as undirected.

 If the #ACTOR ATTRIBUTES, #VERTEX ATTRIBUTES or #EDGE ATTRIBUTES sections are empty,
 no attributes are created.

 The #LAYERS, #ACTORS and #VERTICES sections are useful only if attributes are present, or if there
 are actors that are not present in any layer (#ACTORS), or if there are isolated vertices
 (#VERTICES), otherwise they can be omitted.

 If no section is specified, #EDGES is the default.

 Therefore, a non attributed, undirected multiplex network file can be as simple as:
 \verb{
 ----------------------------------
 Actor1,Actor2,Layer1
 Actor1,Actor3,Layer1
 Actor4,Actor2,Layer2
 ----------------------------------
 }

If interlayer edges exist, then type "multilayer" must be specified, and layers and edges
are formatted in a different way:

\verb{
#VERSION
3.0

#TYPE
multilayer

#ACTOR ATTRIBUTES
AttributeName1,STRING
AttributeName2,NUMERIC
-- etc.

#NODE ATTRIBUTES
LayerName1,AttributeName1,STRING
LayerName1,AttributeName2,NUMERIC
LayerName2,AttributeName3,STRING
-- etc.

#EDGE ATTRIBUTES
-- edge attributes can be defined for specific layers:
LayerName1,AttributeName,STRING
LayerName1,AttributeName,NUMERIC
-- or for all layers (called global attributes):
AnotherAttributeName,NUMERIC
-- etc.

#LAYERS
-- LayerName1,LayerName1,UNDIRECTED
-- LayerName2,LayerName2,DIRECTED
-- LayerName3,LayerName3,DIRECTED,LOOPS
-- LayerName1,LayerName2,DIRECTED
-- etc.
-- all intra-layer specifications (where the first and second layers are the same)
-- should be listed first.
-- LOOPS is only allowed for intra-layer specifications.

#ACTORS
ActorName1,AttributeValueList...
ActorName2,AttributeValueList...
-- etc.

#VERTICES
ActorName1,LayerName1,AttributeValueList...
ActorName1,LayerName2,AttributeValueList...
-- etc.

#EDGES
-- ActorName1,LayerName1,ActorName2,LayerName2,LocalAttributeValueList,GlobalAttributeValueList...
-- etc.

----------------------------------
}
}
\usage{
read_ml(file, name = "unnamed", sep = ',', aligned = FALSE)
write_ml(n, file, format = "multilayer", layers = character(0),
  sep = ',', merge.actors = TRUE, all.actors = FALSE)
}
\arguments{
\item{file}{The path of the file storing the multilayer network.}
\item{name}{The name of the multilayer network.}
\item{n}{A multilayer network.}
\item{layers}{If specific layers are passed to the function, only those layers are saved to file.}
\item{format}{Either "multilayer", to use the package's internal format, or "graphml".}
\item{sep}{The character used in the file to separate text fields.}
\item{aligned}{If \code{true}, all actors are added to all layers.}
\item{merge.actors}{Whether the nodes corresponding to each single actor should be merged into a single node (\code{true}) or kept separated (\code{false}), when \code{format = "graphml"} is used.}
\item{all.actors}{Whether all actors in the multilayer network should be included in the output file (true) or only those present in at least one of the input layers (false), when \code{format = "graphml"} and \code{merge.actors = TRUE} are used.}
}
\value{
\code{read_ml} returns a multilayer network. \code{write_ml} does not return any value.
}
\seealso{
\link{multinet.predefined}, \link{multinet.generation}
}
\examples{
# writing a network to file...
file <- tempfile("aucs.mpx")
net <- ml_aucs()
write_ml(net,file)
# ...and reading it back into a variable
net <- read_ml(file,"AUCS")
net
# the following network has more nodes, because all
# actors are replicated to all graphs
net_aligned <- read_ml(file,"AUCS",aligned=TRUE)
net_aligned
}
