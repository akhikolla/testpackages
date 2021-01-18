\name{plot.simgc}
\alias{plot.simgc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Geostatistical Data Simulated From Gaussian Copula
}

\usage{
\method{plot}{simgc}(x, index = 1, plottype = "Text", xlab = "xloc", ylab = "yloc",
     xlim = NULL, ylim = NULL, pch = 20, textcex = 0.8, plotcex = 1,
     angle = 60, col = 4, col.regions = gray(90:0/100),...)
          }

\description{
Three plots will be generated. A level plot with the number of counts at given locations; a
level plot with point referenced locations and varying colors and a 3-D scatter plot.
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{simgc}, typically generated from function \code{\link{simgc}}.}

  \item{index}{the index of the simulated data, need to be specified since \code{simgc} can
  simulate multiple datasets simultaneously.}

  \item{plottype}{the type of the printed plot, can be "Text", "Dot", or "3D". When plottype =
  "Text", a 2-D plot is generated with exact counts at observed locations. When plottype =
  "Dot", a 2-D dot plot is generated and when plottype = "3D" a 3-D scatter plot is printed.}

  \item{xlab, ylab}{a title for the x and y axis.}

  \item{xlim, ylim}{numeric vectors of length 2, giving the x and y coordinates ranges.
  if they equal to \code{NULL} then they will be adjusted from the data.}

  \item{pch}{plotting character, i.e., symbol to use in the 3-D scatter plot.}

  \item{textcex}{a numerical value giving the amount by which plotting text should be magnified relative to the default.}

    \item{plotcex}{a numerical value giving the amount by which plotting symbols should be magnified relative to the default.}

  \item{angle}{angle between x and y axis.}

  \item{col}{color of the text.}

  \item{col.regions}{color vector to be used reflecting magnitude of the dataset at different locations.
  The general idea is that this should be a color vector of gradually varying color.}

  \item{...}{further arguments passed to plot and panel settings.}
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.mlegc}},
\code{\link{plot.predgc}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Data Visualization}
\keyword{Simulation}
