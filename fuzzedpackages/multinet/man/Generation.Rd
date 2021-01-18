\name{multinet.generation}
\alias{multinet.generation}
\alias{grow_ml}
\alias{evolution_pa_ml}
\alias{evolution_er_ml}
\title{
Generation of multilayer networks
}
\description{
The grow_ml function generates a multilayer network by letting it grow for a number of steps, where for each step three events can happen: (1) evolution according to internal dynamics (in which case a specific internal evolution model is used), (2) evolution importing edges from another layer, and (3) no action.

The functions evolution_pa_ml and evolution_er_ml define, respectively, an evolutionary model based on preferential attachment and an evolutionary model where edges are created by choosing random end points, as in the ER random graph model.
}
\usage{
grow_ml(num.actors, num.steps, models, pr.internal, pr.external, dependency)
evolution_pa_ml(m0,m)
evolution_er_ml(n)
}
\arguments{
\item{num.actors}{The number of actors from which new nodes are selected during the generation process.}
\item{num.steps}{Number of timestamps.}
\item{models}{A vector containing one evolutionary model for each layer to be generated. Evolutionary models are defined using the \code{evolution_*_ml} functions.}
\item{pr.internal}{A vector with (for each layer) the probability that at each step the layer evolves according to the internal evolutionary model.}
\item{pr.external}{A vector with (for each layer) the probability that at each step the layer evolves importing edges from another layer.}
\item{dependency}{A matrix LxL where element (i,j) indicates the probability that layer i will import an edge from layer j in case an external event is triggered.}
\item{m0}{Initial number of nodes.}
\item{m}{Number of edges created for each new vertex joining the network.}
\item{n}{Number of vertices (created at the beginning, before starting adding edges).}
}
\value{
\code{grow_ml} returns a multilayer network. \code{evolution_*_ml} return evolutionary models that are used by \code{grow_ml} to decide how each layer should grow.
}
\seealso{
\link{multinet.predefined}, \link{multinet.IO}
}
\references{
Magnani, Matteo, and Luca Rossi. 2013. Formation of Multiple Networks. In Social Computing, Behavioral-Cultural Modeling and Prediction, 257-264. Springer Berlin Heidelberg.
}
\examples{
# we generate a network with two layers, one growing according
# to the Preferential Attachment model and one growing by selecting
# new edges uniformly at random.
models <- c(evolution_pa_ml(3,1), evolution_er_ml(50))
# all the probability vectors must have the same number of
# fields, one for each layer: two in this example
# by defining pr.internal and pr.external, we are also implicitely defining
# pr.no.action (1 minus the other probabilities, for each field/layer).
pr_external <- c(.5,0)
pr_internal <- c(.5,.5)
# each layer will import edges from the other if needed
# (not the second layer in this example: it has 0 probability of external events)
dependency <- matrix(c(0,1,1,0),2,2)
# 100 steps of network growing, adding actors from a pool of 100
grow_ml(100, 100, models, pr_internal, pr_external, dependency)
}
