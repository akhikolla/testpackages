\name{predict.tfCox}
\alias{predict.tfCox}
\title{Predict for a New Covariate Matrix and fit from \code{\link{tfCox}}
}
\description{This function makes predictions from a specified covariate matrix for a fit of the class "tfCox".
}
\usage{
\method{predict}{tfCox}(object, newX, which.lambda=1, \dots)
}
\arguments{
  \item{object}{
an object of the class "tfCox"
}
  \item{newX}{
a n x p covariate matrix
}
  \item{which.lambda}{
the index for the model of interest to be plotted. \code{which.lambda} corresponds to the model fit in \code{lambda.seq} and should be integer between 1 to \code{length(fit$lambda.seq)}. In other words, the fit from \code{fit$theta.list[[which.lambda]]} will be plotted.
}
  \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}
\details{
Prediction for the new data point is implemented by constant or linear interpolation. 0th order trend filtering will have constant interpolation, and 1th or higher order trend filtering will have linear interpolation.
}
\value{
A n x p matrix containing the fitted theta values.
}

\author{
Jiacheng Wu
}

\seealso{
\code{\link{tfCox}}
}
