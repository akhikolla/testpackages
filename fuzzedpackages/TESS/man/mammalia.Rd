\name{mammalia}
\alias{mammalia}
\title{Dated family level mammalian phylogeny from Meredith et al. (2011): Impacts of the cretaceous terrestrial revolution and kpg extinction on mammal diversification.}
\description{
  This phylogeny describes the species relationship and divergence times of the class Mammalia with 1-3 species included per family, published in Meredith et al. (2011).
}
\usage{
data(mammalia)
}
\format{
  The phylogeny is stored as an object of class \code{"phylo"}.
  The structure is described in the help page of the function
  \code{\link{read.tree}} of the package \code{\link{ape}}.
}
\source{
   Meredith, R. et al. (2011): Impacts of the cretaceous terrestrial revolution and kpg extinction on mammal diversification. Science, 334:521-524
}

\examples{
# load the data
data(mammalia)

# safe the current settings of the plotting device
op <- par()

# set the line width for drawing thinner lines for the branches
par(cex = 0.3)

# plot the mammalian phylogeny
plot(mammalia)

# restore the settings of the device
par(op)
}
\keyword{datasets}
