\name{plotIMGraph}
\alias{plotIMGraph}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plotIMGraph
}
\description{
This function takes a graph object returned from IMaGES (takes the form of a named list containing \code{.graph} and \code{.params}) and plots it. Using plotAll or plotMarkovs is recommended unless you only want to see one specific graph.
}
\usage{
plotIMGraph(graph.object, title="Global")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph.object}{Graph object returned as part of IMaGES run}
  \item{title}{Title of resulting plot of IMGraph}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
  \item{plot }{Prints plot to the screen}
%% ...
}
\references{
  J.D. Ramsey, S.J. Hanson, et al. (2010) Six problems for causal inference
  from fMRI. \emph{NeuroImage} \bold{49}, 1545-1558.
}
\author{
  Noah Frazier-Logue (\email{n.frazier.logue@nyu.edu})
  
  Stephen Jose Hanson (\email{jose@rubic.rutgers.edu})
}

\seealso{
\code{\link{plotMarkovs}}, \code{\link{plotAll}}, \code{\link{IMaGES}}
}
\examples{
## Load predefined data
data(IMData)

#run IMaGES
im.fits <- IMaGES(matrices=IMData, penalty=3, num.markovs=5)

#plot individual graph
plotIMGraph(im.fits$.single.graph[[1]], title="Global")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graph}% use one of  RShowDoc("KEYWORDS")
\keyword{plot}% __ONLY ONE__ keyword per line
