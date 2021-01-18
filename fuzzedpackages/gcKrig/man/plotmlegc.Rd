\name{plot.mlegc}
\alias{plot.mlegc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Geostatistical Data and Fitted Mean
}

\usage{
\method{plot}{mlegc}(x, plotdata = "Observed", plottype = "2D", xlab = "xloc", ylab = "yloc",
     xlim = NULL, ylim = NULL, pch = 20, textcex = 0.8, plotcex = 1,
     angle = 60, col = 4, col.regions = gray(90:0/100),...)
          }

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{mlegc} inherited from function \code{\link{mlegc}}.
}

  \item{plotdata}{the data to be plotted. Can be either "Observed" if the original counts are used or
  "Fitted" if the fitted mean at different locations are used.}

  \item{plottype}{the type of the plot. Can be either "2D" for 2-D level plot or "3D" for 3-D scatterplot.}

  \item{xlab, ylab}{a title for the x and y axis.}

  \item{xlim, ylim}{numeric vectors of length 2, giving the x and y coordinates ranges.
  if they equal to \code{NULL} then they will be adjusted from the data.}

  \item{pch}{plotting character, i.e., the symbol to use in the 3-D scatter plot.}

  \item{textcex}{a numerical value giving the amount by which plotting text should be magnified relative to the default.}

    \item{plotcex}{a numerical value giving the amount by which plotting symbols should be magnified relative to the default.}

  \item{angle}{angle between x and y axis.}

  \item{col}{color of the text.}

  \item{col.regions}{color vector to be used reflecting magnitude of the dataset at different locations.
  The general idea is that this should be a color vector of gradually varying color.}

  \item{...}{further arguments passed to plot and panel settings.}



}
\description{
Four plots can be generated: the 2-D level plot or 3-D scatterplot with the number of counts or fitted values.

}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.simgc}},
\code{\link{plot.predgc}}
}


\keyword{Data Visualization}
\keyword{Estimation}
