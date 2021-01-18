\name{plot.contourLines}
\alias{plot.contourLines}
\title{
Plot contour lines
}
\description{
Plot contour lines from list produced by \code{contourLines} function.
}
\usage{

\method{plot}{contourLines}(x, begin=1, end = length(x), add = FALSE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The list of contour lines (created by \code{contourLines}) you want to plot.
}
  \item{begin}{
Beginning position in list of contour lines you want to plot.
}
  \item{end}{
Ending position in list of contour lines you want to plot.
}
  \item{add}{
A boolean value indicating whether the contour lines should be added to an existing plot (add = TRUE) or should be plotted on a new plot (add = FALSE).
}
  \item{\dots}{
Additional arguments that will be passed to the \code{plot} or \code{lines} function.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
This function does not return anything; it only creates a new plot or modifies an existing plot.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Joshua French
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(volcano)
x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
cL <- contourLines(x, y, volcano)
plot.contourLines(cL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

