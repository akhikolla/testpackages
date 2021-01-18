\name{multinet.communities}
\alias{multinet.communities}
\alias{abacus_ml}
\alias{clique_percolation_ml}
\alias{glouvain_ml}
\alias{glouvain2_ml}
\alias{infomap_ml}
\alias{get_community_list_ml}
\alias{modularity_ml}
\title{
Community detection algorithms
}
\description{
Various algorithms to compute communities in multiplex networks, based on frequent itemset mining (abacus), adjacent cliques (clique percolation), modularity optimization (generalized louvain) and random walks (infomap). glouvain2_ml is a more efficient implementation of the original glouvain_ml, no longer based on matrices. get_community_list_ml is a commodity function translating the result of these algorithms into a list of vertex identifiers, and is internally used by the plotting function.}
\usage{
abacus_ml(n, min.actors=3, min.layers=1)
clique_percolation_ml(n, k=3, m=1)
glouvain_ml(n, gamma=1, omega=1, limit=0)
glouvain2_ml(n, omega=1)
infomap_ml(n, overlapping=FALSE, directed=FALSE, self.links=TRUE)

modularity_ml(n, comm.struct, gamma=1, omega=1)
get_community_list_ml(comm.struct, n)
}
\arguments{
\item{n}{A multilayer network.}
\item{min.actors}{Minimum number of actors to form a community.}
\item{min.layers}{Minimum number of times two actors must be in the same single-layer community to be considered in the same multi-layer community.}
\item{k}{Minimum number of actors in a clique. Must be at least 3.}
\item{m}{Minimum number of common layers in a clique.}
\item{gamma}{Resolution parameter for modularity in the generalized louvain method.}
\item{omega}{Inter-layer weight parameter in the generalized louvain method.}
\item{limit}{Limit parameter in the generalized louvain method: if the number of nodes exceeds this limit, then the modularity is computed on the fly without keeping the full data in memory.}
\item{overlapping}{Specifies if overlapping clusters can be returned.}
\item{directed}{Specifies whether the edges should be considered as directed.}
\item{self.links}{Specifies whether self links should be considered or not.}
\item{comm.struct}{The result of a community detection method.}
}
\value{
All community detection algorithms return a data frame where each row contains actor name, layer name and community identifier.

\code{get_community_list_ml} transforms the output of a community detection function into a list by grouping all the nodes having the same community identifier and the same layer.
}
\references{
\itemize{
\item Berlingerio, Michele, Pinelli, Fabio, and Calabrese, Francesco (2013). ABACUS: frequent pAttern mining-BAsed Community discovery in mUltidimensional networkS. Data Mining and Knowledge Discovery, 27(3), 294-320. (for abacus_ml())
\item Afsarmanesh, Nazanin, and Magnani, Matteo (2018). Partial and overlapping community detection in multiplex social networks. Social informatics (for clique_percolation_ml())
\item Mucha, Peter J., Richardson, Thomas, Macon, Kevin, Porter, Mason A., and Onnela, Jukka-Pekka (2010). Community structure in time-dependent, multiscale, and multiplex networks. Science (New York, N.Y.), 328(5980), 876-8. Data Analysis, Statistics and Probability; Physics and Society. (for glouvain_ml())
\item De Domenico, M., Lancichinetti, A., Arenas, A., and Rosvall, M. (2015)
Identifying Modular Flows on Multilayer Networks Reveals Highly Overlapping Organization in Interconnected Systems. PHYSICAL REVIEW X 5, 011027 (for infomap_ml())
}
}
\seealso{\link{multinet.plotting}}
\examples{
net <- ml_florentine()
abacus_ml(net)
clique_percolation_ml(net)
glouvain_ml(net)
glouvain2_ml(net)
infomap_ml(net)
}
