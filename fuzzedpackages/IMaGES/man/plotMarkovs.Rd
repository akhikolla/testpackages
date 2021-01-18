\name{plotMarkovs}
\alias{plotMarkovs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plotMarkovs
}
\description{
This function takes the object returned by an IMaGES run and plots the global structure with its SEM (structural equation modeling) data, as well as the SEM data for each Markov Equivalence Class (MEC) and their respective structures. The function determines the dimensions that most closely represent a square and plots the graphs in that fashion.
}
\usage{
## Load predefined data
data(IMData)

plotMarkovs(im.fits)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{im.fits}{Object returned by IMaGES run
}
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
\code{\link{plotIMGraph}}, \code{\link{plotAll}}, \code{\link{IMaGES}}
}
\examples{
## Load predefined data
data(IMData)

#run IMaGES
im.fits <- IMaGES(matrices=IMData, penalty=3, num.markovs=5)

#plot global graph and all individual graphs with own SEM data
plotMarkovs(im.fits)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graph}% use one of  RShowDoc("KEYWORDS")
\keyword{plot}% __ONLY ONE__ keyword per line
