\name{plotp}
\alias{plotp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot estimates for the parameter that accounts for possible missing data
}
\description{
Plotting command for use on an object of class "indelmiss".
}
\usage{
plotp(x, model, ci = TRUE, cil = 95, ...)
}
\arguments{
 \item{x}{An object of class "indelmiss".}
  \item{model}{One of "M2" or "M4". }
  \item{ci}{TRUE plots confidence intervals around the estimates.}
  \item{cil}{Confidence interval level.}
  \item{\dots}{Any further commands to plot.}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@binghamton.edu}}
\seealso{
See also \code{\link{plot.indelmiss}}.
}
\examples{
indel <- indelrates(datasource = "simulation", seed = 1, taxa = 5,
             mu = 1, nu = 5, phyl = 5000,  pmiss = 0, toi = 1, 
             bgtype="ancestornodes", bg = c(7, 9), 
             zerocorrection = TRUE, rootprob = "stationary", 
             modelnames = c("M1", "M2", "M3", "M4"),
             optmethod = "nlminb")
print(indel)
plotp(indel, model="M4") 
}