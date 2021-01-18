\name{print.indelmiss}
\alias{print.indelmiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print summary information from fit
}
\description{
Summary command for use on an object of class "indelmiss". Depending on the model, the rates (mu: deletion; nu: insertion), missing data proportion (p), and prior probability of gene family absence at the root are printed. If branch groupings (or clades) were specified, then the rates (and corresponding standard errors) are displayed in a matrix with the columns representing the different branch groupings (ordered by the subsets of x$bg where x is an object of class "indelmiss"). The rows represent the gene deletion and insertion rate, respectively.
}
\usage{
\method{print}{indelmiss}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x    }{An object of class "indelmiss".}
  \item{\dots}{Ignore this.}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@binghamton.edu}}
\seealso{
See also \code{\link{indelrates}} and \code{\link{plot.indelmiss}}.
}
\examples{
indel <- indelrates(datasource = "simulation", seed = 1, taxa = 5,
             mu = 1, nu = 3, phyl = 5000,  pmiss = c(0.2, 0.5), toi = c(1, 3), 
             zerocorrection = TRUE, 
             modelnames = c("M1", "M2", "M3", "M4"),
             optmethod = "nlminb")
print(indel)
}