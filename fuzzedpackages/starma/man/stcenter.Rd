\name{stcenter}
\alias{stcenter}
\title{Space-time centering and scaling function}
\description{
\code{stcenter} centers and scales the space-time series \code{data} such that its mean is 0 and its standard error 1.
}
\usage{
stcenter(data, center=TRUE, scale=TRUE)
}
\arguments{
\item{data}{ a matrix or data frame containing the space-time series: row-wise should be the temporal observations, with each column corresponding to a site. }
\item{center}{ a logical value indicating whether the series should be centered or not (subtracting the mean). }
\item{scale}{ a logical value indicating whether the series should be scaled or not (dividing by the empiric stand deviation). }
}
\details{
To be able to apply the three-stage iterative model building procedure method for STARMA models, data must be centered beforehand (since \code{\link{starma}} doesn't estimate an intercept coefficient).

The only difference with the R function \code{\link{scale}} is that it doesn't center and scale column by column, but globally, since all the observations come from a single process in the case of space time series.
}
\value{
An object of the same class as \code{data}, that is either a \code{matrix} or a \code{data.frame}.
}
\author{Felix Cheysson}
\examples{
data <- matrix(rnorm(9400, mean=5, sd=2), 100, 94)
data <- stcenter(data)	# Center and scale the dataset

# Check for mean
sum(data) / (nrow(data) * ncol(data))

# Check for sd
sqrt( sum(data^2) / (nrow(data) * ncol(data) - 1) )
}
\keyword{starma}
\keyword{stcenter}
\keyword{stscale}
\keyword{scale}
