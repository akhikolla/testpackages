\name{multinet.plotting}
\alias{multinet.plotting}
\alias{plot.Rcpp_RMLNetwork}
\alias{plot.ml}
\alias{plot.multinet}
\alias{values2graphics}
\title{
Drawing a multilayer network
}
\description{
The plot function draws a multilayer network. values2graphics is a support function translating discrete attribute values to graphical parameters. values2graphics returns an object with fields legend.text, legend.pch, legend.col and color or shape.
}
\usage{
\S3method{plot}{Rcpp_RMLNetwork}(x,
layout = NULL, grid = NULL, mai = c(.1,.1,.1,.1),
layers = NULL,
vertex.shape = 21, vertex.cex = 1, vertex.size = vertex.cex, vertex.color = 1,
vertex.labels = NULL, vertex.labels.pos = 3,
vertex.labels.offset = .5, vertex.labels.cex = 1, vertex.labels.col=1,
edge.type = 1, edge.width = 1, edge.col = 1, edge.alpha=.5,
edge.arrow.length = 0.1, edge.arrow.angle = 20,
legend.x = NULL, legend.y = NULL,
legend.pch = 20, legend.cex = 0.5,
legend.inset = c(0, 0),
com = NULL, com.cex = 1,
show.layer.names=TRUE, layer.names.cex=1, ...)

values2graphics(values, output = "color")
}
\arguments{
\item{x}{A multilayer network.}
\item{layout}{A data frame indicating the position of nodes. If NULL, the function layout.multiforce.ml is used to compute it.}
\item{grid}{A vector of size 2 indicating the number of rows and columns where to draw the layers.}
\item{mai}{Percentage of each frame reserved as internal margin (left, top, right, bottom). This only concerns vertices: text labels can be printed inside the margin or even outside the frame depending on their offset.}

\item{layers}{A vector of layer names, that determine which layers and in which order are plotted.}

\item{vertex.shape}{Symbol to use for nodes, corresponding to the parameter {pch} of the R \code{points} function.  This can either
be a single character or an integer code for one of a set of
graphics symbols.  See \code{?points} for more details.}
\item{vertex.size}{synonim of \code{vertex.cex}.}
\item{vertex.cex}{Numeric *c*haracter *ex*pansion factor; multiplied by
\code{par("cex")} yields the final node size.}
\item{vertex.color}{Color of the vertexes. If \code{NULL}, all vertexes in the same layer are plotted using the same color.}

\item{vertex.labels}{A character vector or expression specifying the text to be
written besides each node. It corresponds to the parameter \code{labels} of the R \code{text} function.}
\item{vertex.labels.pos}{A position specifier for the text. Values of '1', '2', '3' and
'4', respectively indicate positions below, to the left of,
above and to the right of the specified coordinates. It corresponds to the parameter \code{pos} of the R \code{text} function.}
\item{vertex.labels.offset}{When \code{vertex.labels.pos} is specified, this value gives the offset of the label from the specified coordinate in fractions of a
character width. It corresponds to the parameter \code{offset} of the R \code{text} function.}
\item{vertex.labels.cex}{Numeric *c*haracter *ex*pansion factor; multiplied by
'par("cex")' yields the final character size.  'NULL' and
'NA' are equivalent to '1.0'. It corresponds to the parameter \code{cex} of the R \code{text} function.}
\item{vertex.labels.col}{Color of the labels.}

\item{edge.type}{Edge line type, corresponding to the 'lty' parameter of the R par function.  Line types can either be specified as an
integer (0=blank, 1=solid (default), 2=dashed, 3=dotted,
4=dotdash, 5=longdash, 6=twodash) or as one of the character
strings '"blank"', '"solid"', '"dashed"', '"dotted"',
'"dotdash"', '"longdash"', or '"twodash"', where '"blank"'
uses 'invisible lines' (i.e., does not draw them). See \code{?par} for more details. It accepts a vector of values which are recycled.}
\item{edge.width}{Edge line width, corresponding to the 'lwd' parameter of the R 'par' function. See \code{?par} for more details. It accepts a vector of
values which are recycled.}
\item{edge.col}{Color of the edges.}
\item{edge.alpha}{Transparency of the edges.}

\item{edge.arrow.length}{Length of the edges of the arrow head (in inches) - corresponding to the parameter of the R \code{arrows} function with the same name.}
\item{edge.arrow.angle}{Angle from the shaft of the arrow to the edge of the arrow
head - corresponding to the parameter of the R \code{arrows} function with the same name.}

\item{legend.x, legend.y}{the x and y co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by \link{xy.coords}.}
\item{legend.pch}{the plotting symbols appearing in the legend, as numeric vector or a vector of 1-character strings}
\item{legend.cex}{character expansion factor relative to current par("cex"). Used for text.}
\item{legend.inset}{inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword.}

\item{com}{The result of a community detection algorithm. When this parameter is set, a colored area is added behind each community.}
\item{com.cex}{Increases (>1) or decreases (<1) the margin around the nodes when the colored areas are drawn around the communities.}


\item{show.layer.names}{if TRUE adds the name of each layer at the center bottom of it.}
\item{layer.names.cex}{Increases (>1) or decreases (<1) the size of the layer names.}

\item{...}{Other graphical parameters.}


\item{values}{A vector with values.}

\item{output}{The type of graphical objects the values should be translated to. It can currently take values "color" and "shape".}
}
\seealso{
\link{multinet.layout}, \link{multinet.communities}
}
\examples{
net <- ml_florentine()
\dontrun{
plot(net)
c <- clique_percolation_ml(net)
plot(net, vertex.labels.cex=.5, com=c)
}
net <- ml_aucs()
\dontrun{
plot(net, vertex.labels=NA)
title("AUCS network")
}
values2graphics(c("a", "b", "b", "c"))
}
