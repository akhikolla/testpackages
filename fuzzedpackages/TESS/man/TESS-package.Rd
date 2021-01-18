\name{TESS-package}
\alias{TESS-package}
\alias{TESS}
\docType{package}
\title{Diversification rate estimation and fast simulation of reconstructed phylogenetic trees under tree-wide time-heterogeneous birth-death processes including mass-extinction events}
\description{Simulation of reconstructed phylogenetic trees under tree-wide time-heterogeneous birth-death processes and estimation of parameters under the same model. Speciation and extinction rates can be any function of time and mass-extinction events at specific times can be provided. Trees can be simulated either conditioned on the number of species, the time of the process, or both. Additionally, the likelihood equations are implemented for convenience and can be used for Maximum Likelihood (ML) estimation and Bayesian inference.}
\details{
\tabular{ll}{
Package: \tab TESS\cr
Type: \tab Package\cr
Version: \tab 2.1.0\cr
Date: \tab 2015-10-23\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
}
\author{Sebastian Hoehna and Michael R. May

Maintainer: Sebastian Hoehna <Sebastian.Hoehna@gmail.com>
}
\references{
S. Hoehna: Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374.

S. Hoehna: Likelihood inference of non-constant diversification rates with incomplete taxon sampling. 2014, PLoS One, 9(1), e84184.

S. Hoehna: The time-dependent reconstructed evolutionary process with a key-role for mass-extinction events. 2015, Journal of Theoretical Biology, 380, 321-331.

S. Hoehna, MR May and BR Moore: TESS: Bayesian inference of lineage diversification rates from (incompletely sampled) molecular phylogenies in R. 2015, Bioinformatics.

MR May, S. Hoehna, and BR Moore: A Bayesian approach for detecting mass-extinction events when rates of lineage diversification vary. 2015, Systematic Biology
}
%\keyword{tree, phylogeny, simulation}
\seealso{
\code{ape}
\code{coda}
%\code{\link[ape]{ape}}
}
