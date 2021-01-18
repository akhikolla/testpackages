\name{plot.indelmiss}
\alias{plot.indelmiss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot parameter estimates from the model fit
}
\description{
Plotting command for use on an object of class "indelmiss". Will draw two plots: one with the estimates for the rates and the other for the "missingness" parameter. plot.indelmiss() calls plotp and plotrates.
}
\usage{
\method{plot}{indelmiss}(x, model = NULL, ci = TRUE, cil = 95, ...)
}
\arguments{
  \item{x}{An object of class "indelmiss".}
  \item{model}{One of "M1", "M2", "M3", or "M4". }
  \item{ci}{TRUE plots confidence intervals around the estimates.}
  \item{cil}{Confidence interval level.}
  \item{\dots}{Any further commands to plot.}
}
\author{
Utkarsh J. Dang and G. Brian Golding

\email{udang@binghamton.edu}}
\seealso{
See also \code{\link{indelrates}}, \code{\link{plotrates}} and \code{\link{plotp}}.
}
\examples{
indel <- indelrates(datasource = "simulation", seed = 1, taxa = 5,
             mu = 1, nu = 5, phyl = 5000,  pmiss = 0, toi = 1, 
             bgtype="ancestornodes", bg = c(7, 9), 
             zerocorrection = TRUE, rootprob="maxlik",  
             modelnames = c("M1", "M2", "M3", "M4"),
             optmethod = "nlminb")
print(indel)
plot(indel, model="M4") 
}