\name{plot.sfa}
\alias{plot.sfa}
\title{
Plotting output from sparse factor analysis.
}
\description{
Function for plotting coefficients from sparse factor analysis.  Returns one of two types of plots: either a summary of estimated dimensionality or a scatterplot of two dimensions.
}
\usage{
\method{plot}{sfa}(x,...)
}
\arguments{
  \item{x}{
Object from output of class \env{sfa}.
}
\item{...}{Additional items to pass to \env{plot}.  Options below.}
}
\details{
The function  produces either a summary of estimated dimensionality or a scatter plot of points between two dimensions.  See options below.

\env{type}{ Either \env{"dim"} or \env{"scatter"}.  Whether to return a summary of dimensionality or scatterplot between two dimensions.
}

\env{main, ylabel, xlabel}{ Main title, y axis label, and x axis label.}

\env{dims.scatter}{ Which two dimensions to plot against each other.}

\env{scatter.by}{ Either \env{row} or \env{col}.  Whether to produce scatterplot of two rows or of two columns.}

\env{topbottom}{ Number of points to label at the extreme of each dimension.  Only used for 
\env{type="scatter"}.}
}
\author{
Marc Ratkovic 
}
\references{
In Song Kim, John Londregan, and Marc Ratkovic.  2015.  "Voting, Speechmaking, and the Dimensions of Conflict in the US Senate."  Working paper.
}

\seealso{\link{sfa}}
\examples{
##  See sfa() for example.
}
