\name{feet}
\alias{feet}
\docType{data}
\title{Foot Measurements for Fourth Grade Children}
\description{This data set was collected by the first author in a fourth grade classroom in Ann Arbor, MI, October 1997. We use the shapereg function to make a shape-restricted fit to this data set. "Width" is a continuous response variable, "length" is a continuous predictor variable, and "sex" is a categorical covariate. The constraint is that "width" is increasing with respect to "length".}
\usage{data(feet)}
\format{
  A data frame with 39 observations on the following 8 variables.
  \describe{
    \item{\code{name}}{First name of child.}
    \item{\code{month}}{Birth month.}
    \item{\code{year}}{Birth year.}
    \item{\code{length}}{Length of longer foot (cm).}
    \item{\code{width}}{Width of longer foot (cm), measured at widest part of foot.}
    \item{\code{sex}}{Boy or girl.}
    \item{\code{foot}}{Foot measured (right or left).}
    \item{\code{hand}}{Right- or left-handedness.}
  }
}
\source{
Meyer, M. C. (2006) Wider Shoes for Wider Feet? \emph{Journal of Statistics Education \bold{Volume 14, Number 1}}. 
}
\examples{
    data(feet)
    l <- feet$length
    w <- feet$width
    s <- feet$sex
    plot(l, w, type = "n", xlab = "Foot Length (cm)", ylab = "Foot Width (cm)")
    points(l[s == "G"], w[s == "G"], pch = 24, col = 2)
    points(l[s == "B"], w[s == "B"], pch = 21, col = 4)
    legend("topleft", bty = "n", c("Girl", "Boy"), pch = c(24, 21), col = c(2, 4))
    title("Kidsfeet Width vs Length Scatterplot")
}

\keyword{data sets}
