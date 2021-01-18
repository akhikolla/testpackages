\name{plotgc}
\alias{plotgc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Geostatistical Count Data
}
\description{
This funtion generates two plots describing a geostatistical count data.
The first plot is a bubble plot with size proportional to the response.
The second plot is a lattice plot with text
describing the number of counts.
}
\usage{
plotgc(data = NULL, locs = NULL, bdry = NULL, col = 2, pch = 1,
       textcex = 1, col.regions = gray(90:0/100), size=c(0.3, 2.7), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{the geostatistical count response.}
  \item{locs}{a n by 2 matrix or data frame that indicates the coordinates of locations.}
  \item{bdry}{a list containing the coordinates of boundaries.}
  \item{col}{the color used for response variable in both plots.}
  \item{pch}{the shape used for response variable in the first plot.}
  \item{textcex}{a numerical value giving the amount by which plotting text should be magnified relative to the default.}
  \item{col.regions}{color vector to be used reflecting magnitude of the dataset at different locations.
  The general idea is that this should be a color vector of gradually varying color.}

  \item{size}{the minimum and maximum of the sizes in the first plot. }
  \item{...}{other parameters that control the plotting.}

%%     ~~Describe \code{x} here~~
}



\author{
Zifei Han \email{hanzifei1@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.simgc}}
\code{\link{plot.mlegc}},
\code{\link{plot.predgc}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data Visualization}

