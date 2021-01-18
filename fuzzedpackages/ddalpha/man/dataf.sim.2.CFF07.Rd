\name{dataf.sim.2.CFF07}
\alias{dataf.sim.2.CFF07}
\title{
Model 2 from Cuevas et al. (2007)
}
\description{
Model 2 from Cuevas et al. (2007)

Processes: \cr
X(t) = m_0(t) + e(t), m_0(t) = 30*(1-t)*t^2 + 0.5*abs(sin(20*pi*t)) \cr
Y(t) = an 8-knot spline approximation of X \cr
e(t): Gaussian with mean = 0, cov(X(s), X(t)) = 0.2*exp(-abs(s - t)/0.3)\cr
the processes are discretized at \code{numDiscrets} equally distant points on [0, 1].
}
\usage{
dataf.sim.2.CFF07(numTrain = 100, numTest = 50, numDiscrets = 51, plot = FALSE)
}

\arguments{
  \item{numTrain}{
number of objects in the training sample
  }
  \item{numTest}{
number of objects in the test sample
  }
  \item{numDiscrets}{
number of points for each object
  }
  \item{plot}{
if TRUE the training sample is plotted
  }
}

\format{
  A data strusture containing \code{$learn} and \code{$test} functional data.
  The functional data are given as data structures.
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates.
      \describe{
        \item{\code{args}}{a numeric vector}
        \item{\code{vals}}{a numeric vector}
       }
    }
    \item{\code{labels}}{The classes of the objects: 0 for X(t), 1 for Y(t)}
  }
}

\source{

Cuevas, A., Febrero, M. and Fraiman, R. (2007). Robust estimation and classification for functional data via projection-based depth notions. Computational Statistics 22 481-496.

}
\seealso{
\code{\link{dataf.*}} for other functional data sets

\code{\link{plot.functional}} for building plots of functional data
}
\examples{
## load the dataset
dataf = dataf.sim.2.CFF07(numTrain = 100, numTest = 50, numDiscrets = 51)
learn = dataf$learn
test = dataf$test

## view the classes
unique(learn$labels)

## access the 5th point of the 2nd object
learn$dataf[[2]]$args[5]
learn$dataf[[2]]$vals[5]

\dontrun{
## plot the data
plot(learn)
plot(test)

## or
dataf = dataf.sim.2.CFF07(numTrain = 100, numTest = 50, numDiscrets = 51, plot = TRUE)
}

}
\keyword{datasets}
\keyword{functional}
