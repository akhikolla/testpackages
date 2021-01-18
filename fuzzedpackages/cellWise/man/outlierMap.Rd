\name{outlierMap}
\alias{outlierMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the outlier map.
}
\description{
The outlier map is a diagnostic plot for the output of \code{\link{MacroPCA}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
outlierMap(res,title="Robust PCA",col="black", pch=16,labelOut=TRUE,id=3,
xlim = NULL, ylim = NULL, cex = 1, cex.main=1.2, cex.lab=NULL, cex.axis=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
A list containing the orthogonal distances (\code{OD}), the score distances (\code{SD}) and their respective cut-offs (\code{cutoffOD} and \code{cutoffSD}). Can be the output of \code{\link{MacroPCA}},
rospca::robpca, rospca::rospca.
}
  \item{title}{
  Title of the plot, default is "Robust PCA".
}
  \item{col}{
Colour of the points in the plot, this can be a single colour for all points or a vector or list specifying the colour for each point. The default is "black".
}
  \item{pch}{
Plotting characters or symbol used in the plot, see points for more details. The default is 16 which corresponds to filled circles.
}
  \item{labelOut}{
Logical indicating if outliers should be labelled on the plot, default is \code{TRUE}.
}
  \item{id}{
Number of OD outliers and number of SD outliers to label on the plot, default is 3.
}
\item{xlim}{
Optional argument to set the limits of the \code{x}-axis.
}
\item{ylim}{
Optional argument to set the limits of the \code{y}-axis.
}
\item{cex}{
Optional argument determining the size of the plotted points. See \code{\link[graphics]{plot.default}} for details.
}
\item{cex.main}{
Optional argument determining the size of the main title. See \code{\link[graphics]{plot.default}} for details.
}
\item{cex.lab}{
Optional argument determining the size of the labels. See \code{\link[graphics]{plot.default}} for details.
}
\item{cex.axis}{
Optional argument determining the size of the axes. See \code{\link[graphics]{plot.default}} for details.
}
}
\details{
The outlier map contains the score distances on the x-axis and the orthogonal distances on the y-axis. To detect outliers, cut-offs for both distances are shown, see Hubert et al. (2005).
}

\references{
Hubert, M., Rousseeuw, P. J., and Vanden Branden, K. (2005). ROBPCA: A New Approach to Robust Principal Component Analysis. \emph{Technometrics}, \bold{47}, 64-79.
}
\author{
P.J. Rousseeuw
}

\seealso{
\code{\link{MacroPCA}}
}
\examples{
# empty for now
}