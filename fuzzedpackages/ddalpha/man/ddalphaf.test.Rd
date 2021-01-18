\name{ddalphaf.test}
\alias{ddalphaf.test}
\title{
Test Functional DD-Classifier
}
\description{
Trains functional DD-classifier on the learning sequence of the data and tests it on the testing sequence.
}
\usage{
ddalphaf.test(learn, learnlabels, test, testlabels, disc.type = c("LS", "comp"), ...)
}
\arguments{
  \item{learn}{
list containing lists (functions) of two vectors of equal length, named "args" and "vals": arguments sorted in ascending order and corresponding them values respectively
}
  \item{learnlabels}{
list of output labels of the functional observations
}
  \item{test}{
the testing sequence. Has the same format as \code{learn}
}
  \item{disc.type}{
type of the used discretization scheme. "LS" for \code{\link{ddalphaf.train}}, "comp" for  for \code{\link{compclassf.train}}
}
  \item{testlabels}{
list of output labels of the functinal observations
}
  \item{\dots}{
additional parameters passed to \code{\link{ddalphaf.train}}
}
}

\value{

  \item{error}{
  the part of incorrectly classified data
  }
  \item{correct}{
  the number of correctly classified objects
  }
  \item{incorrect}{
  the number of incorrectly classified objects
  }
  \item{total}{
  the number of classified objects
  }
  \item{ignored}{
  the number of ignored objects (outside the convex hull of the learning data)
  }
  \item{n}{
  the number of objects in the testing sequence
  }
  \item{time}{
  training time
  }

}


\seealso{
\code{\link{ddalphaf.train}} to train the functional DD\eqn{\alpha}-classifier, 
\code{\link{ddalphaf.classify}} for classification using functonal DD\eqn{\alpha}-classifier, 
\code{\link{ddalphaf.getErrorRateCV}} and \code{\link{ddalphaf.getErrorRatePart}} to get error rate of the functional DD-classifier on particular data.
}
\examples{

# load the fdata
df = dataf.growth()

samp = c(35:70)

ddalphaf.test(learn = df$dataf[-samp], learnlabels = df$labels[-samp], 
              test =  df$dataf[samp],  testlabels =  df$labels[samp], 
              adc.args = list(instance = "avr", 
                              numFcn = 2, 
                              numDer = 2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ benchmark }
