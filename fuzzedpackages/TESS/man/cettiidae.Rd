\name{cettiidae}
\alias{cettiidae}
\title{Cettiidae phylogeny from Alstroem et al. (2011)}
\description{
  This phylogeny describes the species relationship and divergence times of the bird family Cettiidae, published in Alstroem et al. (2011).
}
\usage{
data(cettiidae)
}
\format{
  The phylogeny is stored as an object of class \code{"phylo"}.
  The structure is described in the help page of the function
  \code{\link{read.tree}} of the package \code{\link{ape}}.
}
\source{
   Alstroem, P., Hoehna, S., Gelang, M., Ericson, P.G.P, and Olsson, U. (2011) Non-monophyly and intricate morphological evolution within the avian family Cettiidae revealed by multilocus analysis of a taxonomically densely sampled dataset, BMC Evolutionary Biology, 11:352.
}

\examples{
# load the data
data(cettiidae)

# safe the old plotting settings
op <- par()

# set the new plotting settings
par(cex = 0.3)

# plot the phylogeny
plot(cettiidae)

# restore the plotting settings
par(op)
}
\keyword{datasets}
