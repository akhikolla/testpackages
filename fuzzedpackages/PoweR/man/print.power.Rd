\name{print.power}

\alias{print.power}
\alias{print.power1}

\title{Latex table for power simulations}

\description{Transform the power values given by function \code{\link{powcomp.fast}} 
into a LaTeX code for creating the table of power simulations.}

\usage{
\method{print}{power}(x, digits = 3, latex.output = FALSE, template = 1, ...)
}

\arguments{
  \item{x }{power values given by function \code{\link{powcomp.fast}}.}
  \item{digits }{control the number of decimal places. It can take values from 0 to 3.}
  \item{latex.output }{logical. If \code{TRUE}, we output LateX code for the table of power simulations. 
                      If \code{FALSE}, we output this table in the R Console.}
  \item{template }{integer, template to use for the (LaTeX) printing of
  values. Only \code{template = 1} is defined for the moment.}
  \item{... }{further arguments passed to or from other methods.}
}

\references{
Pierre Lafaye de Micheaux, Viet Anh Tran (2016). PoweR: A
Reproducible Research Tool to Ease Monte Carlo Power Simulation
Studies for Studies for Goodness-of-fit Tests in R. \emph{Journal of Statistical Software}, \bold{69(3)}, 1--42. doi:10.18637/jss.v069.i03 

Puig, P. and Stephens, M. A. (2000), Tests of fit for the Laplace distribution, with applications, 
\emph{Technometrics}, \bold{42}, 417--424.
}
\author{P. Lafaye de Micheaux, V. A. Tran}
\seealso{See \code{\link{print.critvalues}}.}

\examples{
## Regenerate Table 6 from Puig (2000) (page 424)
# Change M = 50000 for more accurate results
M <- 10
law.index <- 1
vectn <- c(10,15,20,35,50,75,100)
level <- c(0.05)
stat.indices <- c(43,44,42,45,46)
law.indices <- c(2,3,4)
alter <- list(stat43 = 3,stat44 = 3,stat42 = 3,stat45 = 3,stat46 = 3)
critval <- many.crit(law.index,stat.indices,M,vectn,level,alter,law.pars = NULL,parstat = NULL)
table6 <- powcomp.fast(law.indices,stat.indices,vectn,M,level,critval = critval,alter,
                       parlaws = NULL,parstats = NULL,nbclus = 1)
print.power(table6,digits=0,latex.output = TRUE)
}

\keyword{print}