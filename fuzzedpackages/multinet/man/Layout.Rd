\name{multinet.layout}
\alias{multinet.layout}
\alias{layout_multiforce_ml}
\alias{layout_circular_ml}
\title{
Layouts
}
\description{
These functions compute xyz coordinates for each node in the network.
}
\usage{
layout_multiforce_ml(n, w_in = 1, w_inter = 1, gravity = 0, iterations = 100)
layout_circular_ml(n)
}
\arguments{
\item{n}{A multilayer network.}
\item{w_in}{An array with weights for intralayer forces, or a single number if weights are the same for all layers. When \code{w_in} is positive, vertices in the corresponding layer will be positioned as if a force was applied to them, repelling vertices that are close to each other and attracting adjacent vertices, all proportional to the specified weight.}
\item{w_inter}{An array with weights for interlayer forces, or a single number if weights are the same for all layers. When \code{w_inter} is positive, vertices in the corresponding layer will be positioned as if a force was applied to them, trying to keep them aligned with the vertices corresponding to the same actors on other layers, proportionally to the specified weight.}
\item{gravity}{An array with weights for gravity forces, or a single number if weights are the same for all layers. This parameter results in the application of a force to the vertices, directed toward the center of the plot. It can be useful when there there are multiple components, so that they do not drift away from each other because of the repulsion force applied to their vertices.}
\item{iterations}{Number of iterations.}
}
\references{
Fatemi, Zahra, Salehi, Mostafa, & Magnani, Matteo (2018). A generalised force-based layout for multiplex sociograms. Social Informatics
}
\seealso{
\link{multinet.plotting}
}
\examples{
net <- ml_florentine()
layout_multiforce_ml(net)
l <- layout_circular_ml(net)
\dontrun{
plot(net,layout=l)}
}
