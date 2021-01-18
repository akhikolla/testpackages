\name{invlogit}
\alias{invlogit}
\title{Inverse logistic functions.}

\description{
Inverse-logit function, transforms continuous values to the range (0, 1)
}

\usage{
invlogit(x)
}

\arguments{
  \item{x}{A vector of continuous values}
}

\details{
        The Inverse-logit function defined as: 
        \eqn{logit^{-1}(x) = e^x/(1+e^x)} transforms continuous values to the range (0, 1),
        which is necessary, since probabilities must be between 0 and 1 and maps 
        from the linear predictor to the probabilities. 
}

\value{
  A vector of estimated probabilities
}

\author{Artemis Toumazi \email{artemis.toumazi@gmail.com}, 
        Moreno Ursino \email{moreno.ursino@inserm.fr},
        Sarah Zohar \email{sarah.zohar@inserm.fr}
        }

\keyword{models}
 
