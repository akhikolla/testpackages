\name{humanity_model}
\alias{humanity.X}
\alias{humanity.Y}
\alias{humanity.Xt}
\alias{humanity.Yt}
\docType{data}
\title{
data from the humanity model
}
\description{
This data set provides the training data and testing data from the 'diplomatic and military operations in a non-warfighting
domain' (DIAMOND) simulator. It porduces the number of casualties during the second day to sixth day after the earthquake and
volcanic eruption in Giarre and Catania. See  (Overstall and Woods (2016)) for details.
}
\usage{
	data(humanity_model)
}
\format{
  Four data frame with observations on the following variables.
  \describe{
    \item{\code{humanity.X}}{A matrix of the training inputs.}
    \item{\code{humanity.Y}}{A matrix of the output of the calsualties from the second to sixth day after the the earthquake and
volcanic eruption for each set of training inputs.}
    \item{\code{humanity.Xt}}{A matrix of the test inputs.}
    \item{\code{humanity.Yt}}{A matrix of the test output of the calsualties.}
  }
}
\references{


A. M. Overstall and D. C. Woods (2016). Multivariate emulation of computer simulators: model selection and diagnostics with application to a humanitarian relief model. Journal of the Royal Statistical Society: Series C (Applied Statistics), 65(4):483-505.

B. Taylor and A. Lane. Development of a novel family of military campaign simulation models. Journal of the Operational Research Society, 55(4):333-339, 2004.

}

\keyword{datasets}
