\name{plottree}
\alias{plottree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the tree used the branches colored according to the different specified branch groupings (or clades) following unique rates.
}
\description{
Plotting command for use on an object of class "indelmiss".
}
\usage{
plottree(x, toilabel = TRUE, colors = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "indelmiss".}
  \item{toilabel}{If this is TRUE, a plus sign is printed next to the taxa of interest for which a missing data proportion was estimated. Note that the taxa labels being referred to can be seen by using ape::tiplabels().}
  \item{colors}{Vector of colours the same length as length(x$bg). Note that these colours are used to colour the different branch groupings associated with unique insertion and/or deletion rates.}
  \item{\dots}{Any further commands to ape::plot.phylo.}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@binghamton.edu}}
\seealso{
See also \code{\link{plot.indelmiss}} and \code{\link[ape]{plot.phylo}}.
}
\examples{
indel <- indelrates(datasource = "simulation", seed = 1, taxa = 5,
             mu = 1, nu = 5, phyl = 5000,  pmiss = 0, toi = 1, 
             bgtype="ancestornodes", bg = 7, 
             zerocorrection = TRUE, 
             modelnames = c("M1", "M2", "M3", "M4"),
             optmethod = "nlminb")
print(indel)
plottree(indel,colors=c("blue","red"))
ape::tiplabels()
}