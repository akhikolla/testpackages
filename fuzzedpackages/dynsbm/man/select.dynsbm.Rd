\name{select.dynsbm}
\alias{select.dynsbm}
\title{
Dynamic stochastic block model estimation for different number of groups.
}
\description{
Estimation of dynamic stockastic block models for different number of groups.
Each model combines a stochastic block model (SBM) for its static part with independent Markov chains for the evolution of the nodes groups through time.  
}
\usage{
select.dynsbm(Y, present=NULL, Qmin, Qmax,
			edge.type=c("binary","discrete","continuous"), K=-1,
			directed=FALSE, self.loop=FALSE,
			nb.cores=1,
			iter.max=20, nstart=25, perturbation.rate=0.2,
                        fixed.param=FALSE, bipartition=NULL,
			plot=TRUE)
}
\arguments{
  \item{Y}{An object of class \code{array} of dimension (T x N x N) containing T adjacency matrices of size (N x N), where N is the number of nodes in the network and T is the number of time points.}
  \item{present}{\code{NULL} or an object of class \code{matrix} of size (N x T)
    containing the presence/absence (coded with 1/0 respectively) of
    each N nodes at each of the T time points. When set to \code{NULL}, this object is deduced from \code{Y} (see the "Details" section). Any node must be present at least once among the time points.}
  \item{Qmin}{Minimum number of groups >1.}
  \item{Qmax}{Maximum number of groups.}
  \item{edge.type}{Type of adjacency matrices. This should be (an unambiguous abbreviation of) one of \code{binary}, \code{discrete} or \code{continuous}. 
  See the "Details" section.}
  \item{K}{Only if \code{edge.type=="discrete"}. Number of non-zero discrete values (i.e. in {1,..,K}).} 
  \item{directed}{If \code{TRUE}, the network is supposed to be directed (and
    therefore \code{Y} is supposed to be asymmetric).}
  \item{self.loop}{If \code{TRUE}, self-loops (edges from one node to the same node) are allowed and taken into acount in the estimation procedure.}   
  \item{nb.cores}{Number of cores to use, i.e. how many child
    processes and how many threads will be run simultaneously during the
    initialization and the estimation steps respectively.}
  \item{iter.max}{Maximal number of algorithm iterations.}
  \item{nstart}{Number of starting points for the iterative estimation algorithm.  See the "Details" section.}
  \item{perturbation.rate}{Rate of perturbation (in [0,1], see \code{nstart}) for the iterative estimation algorithm. This rate is
    the fraction of nodes for which its group is randomly shuffled.}
  \item{fixed.param}{If \code{TRUE}, the model parameters remain fixed and
    constant in time. By default, \code{fixed.param} is automatically
    set to \code{TRUE} in the bipartite case (i.e. \code{bipartition} is not
    \code{NULL}; see the "Details" section).}
  \item{bipartition}{\code{NULL} or a vector of size N specifying a node bipartition
    in the case of bipartite networks (see the "Details" section). Each element of this vector is set to 1 or 2 to specify the node belongs to the
    first or second set of nodes.}
  \item{plot}{Display a plot with the loglikelihood and the ICL criteria if \code{edge.type=="binary"/"continuous"}. See the "Details" section.} 
}
\details{
  This function deals with binary or weighted dynamic/temporal/evolving networks (with discrete or continuous edges).
  The adjacency matrices must be coded with 0/1 in the \code{binary} case, 
  with 0/y where y belongs to the set {1,..,K} in the \code{discrete} case 
  or with 0/y where y is numeric, must be positive and is supposed to fit a gaussian mixture in
  the \code{continuous} case.

  Presence/absence information allows to model node's arrival or
  departure, birth or death, or simply enables to specify missing data
  (as absent nodes). If this information is missing (NULL), the presence/absence is deduced automatically from \code{Y} 
    by searching for nodes that do not participate in any edges (lines/columns of O in \code{Y}) and declaring them as absent. 
  This function does not support the existence of nodes that are never present (error message in this case).  
  
  The estimation algorithm is iterative and rely on a starting point. Therefore, it is possible to start the algorithm many times with 'nstart' starting points. 
  The first starting point is obtained with an ad-hoc use of the \code{kmeans} function. 
  The follwing starting point are obtained by perturbating the first one (see \code{perturbation.rate}).
  The greater \code{nstart}, the more accurate the results.
  
  To select the best number of groups, the "elbow" method consists in finding the point where the slope of the loglikelihood is decreasing (i.e. the loglikelihood is reaching a plateau).
  If \code{edge.type=="binary"}, the ICL criteria (plotted in red) has
  to be used : the best number of nodes is supposed to maximize the ICL
  criteria. 

  This function has been extended to the case of bipartite networks. In
  this case, despite \code{Y} has to be of dimension (T x N x N), it is
  possible to give a bipartition of the nodes into two disjoint
  sets. For statistical reasons, \code{fixed.param} is automatically set
  to \code{TRUE}. Given the total number of groups Q between Qmin and
  Qmax, there is Q/2 groups for each set of nodes (when Q is odd,
  there is floor(Q/2)+1 groups for the largest set of nodes);
  however, there is no guaranty that the final groups are
  coherent with the bipartition, i.e. that any group is composed by
  nodes of one of the two sets (if not, a warning message is generated).
}
\value{
Returns a list of \code{dynsbm} objects. Each object of class \code{dynsbm} is a \code{list} with the following components:
\item{transition}{The Markov chain transition \code{matrix} of size (Q x Q).}
\item{membership}{An object of class \code{matrix} of size (N x T)
  containing the group membership estimated by MAP (>0, =0 for absent nodes).}
\item{beta}{An object of class \code{matrix} of size (T x Q x Q) containing
  the sparsity parameters of the model.}
\item{gamma}{Only if \code{edge.type=="discrete"}. An object of class \code{matrix} of size (T x Q x Q x K) containing
  the model parameters.}
\item{mu/sigma}{Only \code{edge.type=="continuous"}An object of class \code{matrix} of size (T x Q x Q)
  and a vector containing the model parameters.}
\item{loglikelihood}{Completed data log-likelihood.}
\item{iter}{Number of used algorithm iterations.}
\item{directed}{Specifies whether the model is build for directed networks.}
\item{self.loop}{Specifies whether the model allows self-loops.}
}
\references{
Catherine Matias and Vincent Miele, \emph{Statistical clustering of temporal networks through a dynamic stochastic block model}, Journal of the Royal Statistical Society: Series B (2017)
http://dx.doi.org/10.1111/rssb.12200
http://arxiv.org/abs/1506.07464

Vincent Miele and Catherine Matias, \emph{Revealing the hidden structure
  of dynamic ecological networks},  Royal Society Open Science (2017)
http://dx.doi.org/10.1098/rsos.170251
https://arxiv.org/abs/1701.01355
}
\author{
Authors: Catherine Matias, Vincent Miele
  
Maintainer: Vincent Miele <vincent.miele@univ-lyon1.fr>
}
\examples{
data(simdataT5Q4N40binary)

## estimation for Q=1..6 groups
list.dynsbm <- select.dynsbm(simdataT5Q4N40binary, 
				Qmin=1, Qmax=6, edge.type="binary", nstart=1)
				
\dontrun{
## better to use nstart>1 starting points
## but estimation can take 1-2 minutes
list.dynsbm <- select.dynsbm(simdataT5Q4N40binary, 
				Qmin=1, Qmax=6, edge.type="binary", nstart=25)
}

## selection of Q=4
dynsbm <- list.dynsbm[[4]]

## plotting intra/inter connectivity patterns
connectivity.plot(dynsbm, simdataT5Q4N40binary)

## plotting switches between groups
alluvial.plot(dynsbm)
}