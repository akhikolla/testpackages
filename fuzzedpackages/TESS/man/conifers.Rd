\name{conifers}
\alias{conifers}
\title{Conifer phylogeny from Leslie et al. (2012)}
\description{
  This phylogeny describes the species relationships and divergence times of the plant order Conifera, published in Leslie et al. (2012).
}
\usage{
data(conifers)
}
\format{
  The phylogeny is stored as an object of class \code{"phylo"}.
  The structure is described in the help page of the function
  \code{\link{read.tree}} of the package \code{\link{ape}}.
}
\source{
  Leslie, A. B., J. M. Beaulieu, H. S. Rai, P. R. Crane, M. J. Donoghue, and S. Mathews. 2012.  Hemisphere-scale differences in conifer evolutionary dynamics. Proceedings of the National Academy of Sciences 109:16217-16221.
}

\examples{
# load the tree
data(conifers)

# safe the settings of the plotting device
op <- par()

# set the line width for plotting the branches
par(cex = 0.3)

# plot the phylogenetic tree
plot(conifers)

# restore the settings of the device
par(op)
}
\keyword{datasets}
