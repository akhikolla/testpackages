\name{expSBM_init}
\alias{expSBM_init}
\title{expSBM_init}
\description{
Initialization step for the variational expectation maximization algorithm. 
}
\usage{
expSBM_init(edgelist, K, 
   method = c("random", "SBM_binary", "SBM_gaussian", "spectral"), 
   sbm_density = 0.5, blur_value = 1) 
}
\arguments{
  \item{edgelist}{A matrix with 4 columns: on the first column the sender node, on the second the receiver, on the third either a one or zero to indicate whether it is an interaction or a non-interaction respectively, on the fourth the corresponding exponential length.}
  \item{K}{Number of latent groups.}
  \item{method}{Method used to initialise the allocations. Can be one of \code{random}, \code{SBM_binary}, \code{SBM_gaussian} or \code{spectral}. See details.}
  \item{sbm_density}{If \code{method == "SBM_binary"}, this is the target density for the thresholded binary stochastic blockmodel.}
  \item{blur_value}{A value between 0 and 1. If 1, the initialization method returns a hard partition where each node belongs to one group and one only. Reducing this value introduces noise, i.e. it gradually transforms the hard clustering into a soft clustering where each node is equally likely to belong to any of the K given clusters.}
}
\details{
All initialisation methods return a \code{NxK} matrix indicating the partitioning of the nodes.

The method \code{random} intialises the allocation variables uniformly at random.

The method \code{SBM_binary} first calculates the aggregated interaction and non-interaction times for each pair of nodes. Then, it calculates the portion of time when the nodes where interacting over the whole time period. Then it obtains an adjacency matrix by thresholding these values, i.e. values above a given threshold are replaced by ones and values below the threshold are replaced by zeros. The threshold is chosen by setting the parameter \code{sbm_density} which defines the desired density of the graph. Once the adjacency matrix is obtained, a binary stochastic blockmodel is fit on the data hence obtaining the partition.

The method \code{SBM_gaussian} aggregates the interaction values and non-interaction values for each pair of nodes. Then it log-transforms both of these quantities. Then it fits a stochastic blockmodel with multivariate Gaussian edges to obtain the partition.

The method \code{spectral} first calculates the aggregated interaction and non-interaction times for each pair of nodes. Then, it calculates the portion of time when the nodes where interacting over the whole time period. Then it performs model-based clustering on the Laplacian associated to this weighted matrix.
}
\value{
A \code{NxK} matrix indicating the partitioning of the nodes.
}
\examples{
set.seed(12345)
data(high_school)
K <- 4
lambda_init <- rep(1/K,K)
expSBM_init(high_school$edgelist, K, "random")
}


