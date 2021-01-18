\name{saturnin-package}
\alias{saturnin-package}
\alias{saturnin}
\docType{package}
\title{
Spanning Trees Used for Network Inference
}
\description{
Bayesian inference of graphical model structures using spanning trees.  
For further details on the considered framework, we refer
the reader to the paper quoted in the references section.
}
\details{
\tabular{ll}{
Package: \tab saturnin\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-10\cr
License: \tab GPL-2\cr
}
}
\author{
Loïc Schwaller

Maintainer: Loïc Schwaller <loic.schwaller@ens-lyon.org>

}
\references{
This package implements the method described in the paper "Bayesian Inference of Graphical Model Structures Using Trees" by L. Schwaller, S. Robin, M. Stumpf, 2015 (submitted and availavable on arXiv).
}

\keyword{ graphical models }
\keyword{ spanning trees }
\keyword{ hyper Markov }
\keyword{ bayesian inference }
\keyword{ network inference }

\examples{
library('saturnin')
data(data_gaussian)

W <- lweights_gaussian(data_gaussian)
prob <- edge.prob(W, log = TRUE, account.prior = TRUE, q0 = 0.5)
}
