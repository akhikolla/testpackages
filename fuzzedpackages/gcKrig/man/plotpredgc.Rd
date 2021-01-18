\name{plot.predgc}
\alias{plot.predgc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Geostatistical Data at Sampling and Prediction Locations
}



\usage{
\method{plot}{predgc}(x, plottype = "2D", xlab = "xloc", ylab = "yloc", xlim = NULL,
     ylim = NULL, pch = 20, textcex = 0.6, plotcex = 1, angle = 60,
     col = c(2, 4), col.regions = gray(90:0/100),...)
          }

\description{
Five plots can be generated. A level plot with the number of counts at both observed and
prediction locations; a
level plot with predicted means (intensity); a level plot with the predicted counts;
a level plot with estimated variances of the prediction; a
3-D scatter plot with both observed and predicted counts.
}


%- maybe also 'usage' for other objects documented here.
\arguments{

 \item{x}{an object of class \code{predgc} inherited from function \code{\link{predgc}}.}

 \item{plottype}{can be one of the following: "2D", "Predicted Counts", "Predicted Mean",
 "Predicted Variance" or "3D". Default is "2D" which generates a 2-D contour plot with
 both observed and predicted counts. With arguments "Predicted Counts", "Predicted Mean" and
 "Predicted Variance", a 2-D level plot will be generated with the corresponding data.
  When "3D" is used, a 3-D scatter plot will be displayed with observed and predicted counts.
 }

  \item{xlab, ylab}{a title for the x and y axis.}

  \item{xlim, ylim}{numeric vectors of length 2, giving the x and y coordinates ranges.
  if they equal to \code{NULL} then they will be adjusted from the data.}


  \item{pch}{plotting character, i.e., symbol to use in the 3-D scatter plot.}

  \item{textcex}{a numerical value giving the amount by which plotting text should be magnified relative to the default.}

    \item{plotcex}{a numerical value giving the amount by which plotting symbols should be magnified relative to the default.}

  \item{angle}{angle between x and y axis.}

  \item{col}{a numeric vector of length 2 indicating
  color of the plot at sampling and prediction locations.}

  \item{col.regions}{color vector to be used reflecting magnitude of the dataset at different locations.
  The general idea is that this should be a color vector of gradually varying color.}

  \item{...}{further arguments passed to plot and panel settings.}


}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.simgc}},
\code{\link{plot.mlegc}},
\code{\link{mlegc}},
\code{\link{predgc}}
}


\keyword{Data Visualization}
\keyword{Prediction}
