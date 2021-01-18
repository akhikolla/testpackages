\name{SBTMProbs}
\alias{SBTMProbs}
\title{SBTM parameter estimators}
\description{
Evaluates estimates for the model parameters of a Stochastic Block Transition Model. This function can summarise the results obtained through the \code{\link{GreedyICL}} and \code{\link{GreedyMerge}} algorithms.
}

\usage{
SBTMProbs(adj_cube, allocations)
}

\arguments{
  \item{adj_cube}{A binary array of size \code{NxNxT} representing the dynamic network. The generic entry in position \code{[i,j,t]} is equal to \code{1} if \code{i} interacts with \code{j} during the \code{t}-th time frame, or to \code{0} otherwise. Only undirected networks with no self-edges are supported, so each slice of the array must be a symmetric matrix with null elements on the diagonal.}
  \item{allocations}{A matrix of size \code{TxN} denoting the cluster membership of each node at each time. Values should be strictly positive integers. The value zero should be used to identify the inactive nodes.}
}

\value{
  \item{Pi}{A matrix of size \code{(K+1)x(K+1)}. The generic entry \code{[g,h]} denotes the estimated probability that a node moves from group \code{g-1} to group \code{h-1}. The labels are shifted because the allocations take values from \code{0} to \code{K}, with zeros corresponding to inactive nodes. The first row and first column denote the transition probabilities for the group of inactive nodes.}
  \item{Theta}{Connection probability matrix of size \code{KxK} for the SBM part of the model. In this case the group of inactive nodes is omitted. The entry in \code{[g,h]} is the estimated probability that an edge between a node in group \code{g} and a node in group \code{h} is realised. Note that this is only for the case where no information is available regarding the presence or absence of the same edge in the previous time-frame.}
  \item{P}{Edge creation probability matrix of size \code{KxK} for the SBTM part of the model. As in \code{Theta}, the group of inactive nodes is omitted. The entry in \code{[g,h]} is the probability that an edge between a node in group \code{g} and a node in group \code{h} changes its value from \code{0} to \code{1}.}
  \item{Q}{Edge destruction probability matrix of size \code{KxK} for the SBTM part of the model. As in \code{Theta}, the group of inactive nodes is omitted. The entry in \code{[g,h]} is the probability that an edge between a node in group \code{g} and a node in group \code{h} changes its value from \code{1} to \code{0}.}
}

\seealso{
  \code{\link{GreedyICL}}, \code{\link{GreedyMerge}}.
}



