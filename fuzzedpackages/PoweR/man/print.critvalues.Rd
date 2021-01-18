\name{print.critvalues}

\alias{print.critvalues}
\alias{print.critvalues1}

\title{Latex table for critical values}

\description{Transform the critical values given by function \code{\link{many.crit}} 
into a LaTeX code for creating the table of critical values.}

\usage{
\method{print}{critvalues}(x, digits = 3, latex.output =  FALSE, template = 1, ...)
}

\arguments{
  \item{x }{critical values given by function \code{\link{many.crit}}.}
  \item{digits }{integer indicating the number of decimal places to be used.}
  \item{latex.output }{logical. If \code{TRUE}, we output LaTeX code for the table of critical values. 
                      If \code{FALSE}, we output this table in the R Console.}
  \item{template }{integer, template to use for the (LaTeX) printing of
  values. Only \code{template = 1} is defined for the moment.}
  \item{... }{further arguments passed to or from other methods.}
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03


Puig, P. and Stephens, M. A. (2000), Tests of fit for the Laplace distribution, with applications, 
\emph{Technometrics}, \bold{42}, 417--424.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{print.power}}.}

\examples{
## Regenerate Table 1 from Puig (2000) (page 419)
# Take M = 10000 for accurate results
M <- 10
law.index <- 1
vectn <- c(10,15,20,35,50,75,100,1000)
level <- c(0.50,0.25,0.10,0.05,0.025,0.01) 
table1 <- many.crit(law.index,stat.indices = c(43),M,vectn,level,
                    alter = list(stat43=3),law.pars = NULL,parstat = NULL)
print.critvalues(table1,digits=3,latex.output=TRUE)
}

\keyword{print}