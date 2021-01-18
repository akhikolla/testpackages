\name{plot}
\alias{plot.RealVAMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method for RealVAMS
}
\description{
Plot teacher effects and residuals.  The caterpillar plots use a modified version of the plotCI function from R package gplots. According to that package, "Original version [of plotCI] by Bill Venables wvenable@attunga.stats.adelaide.edu.au posted to r-help on Sep. 20, 1997. Enhanced version posted to r-help by Ben Bolker ben@zoo.ufl.edu on Apr. 16, 2001. This version was modified and extended by Gregory R. Warnes greg@warnes.net. Additional changes suggested by Martin Maechler maechler@stat.math.ethz.ch integrated on July 29, 2004."
}
\usage{
\method{plot}{RealVAMS}(x, ..., alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class RealVAMS
}
  \item{\dots}{
other arguments
}
  \item{alpha}{
the significance level for the caterpillar plots
}

}

\value{
Requires user to click window or press "enter" to progress through plots. Returns caterpillar plots (via the package gplots) and residual plots. 
}

\author{
Andrew Karl \email{akarl@asu.edu}, Jennifer Broatch, Jennifer Green, Other authors as listed above for the caterpillar plots.
}


\seealso{
 \code{\link{summary.RealVAMS}}
}
\examples{
data(example.score.data)
data(example.outcome.data)
\donttest{
res<-RealVAMS(example.score.data,example.outcome.data)
plot(res)
}
}
\keyword{ regression }

