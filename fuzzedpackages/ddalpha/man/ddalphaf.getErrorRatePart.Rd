\name{ddalphaf.getErrorRatePart}
\alias{ddalphaf.getErrorRatePart}
\title{
Test Functional DD-Classifier
}
\description{
Performs a benchmark procedure by partitioning the given data. 
On each of \code{times} steps \code{size} observations are removed from the data, the functional DD-classifier is trained on these data and tested on the removed observations.
}
\usage{
ddalphaf.getErrorRatePart(dataf, labels, size = 0.3, times = 10, 
                          disc.type = c("LS", "comp"),  ...)
}
\arguments{
  \item{dataf}{
list containing lists (functions) of two vectors of equal length, named "args" and "vals": arguments sorted in ascending order and corresponding them values respectively
}
  \item{labels}{
list of output labels of the functional observations
}
  \item{size}{
  the excluded sequences size. Either an integer between \eqn{1} and \eqn{n}, or a fraction of data between \eqn{0} and \eqn{1}.
}
  \item{times}{
  the number of times the classifier is trained.
}
  \item{disc.type}{
type of the used discretization scheme. "LS" for \code{\link{ddalphaf.train}}, "comp" for  for \code{\link{compclassf.train}}
}
  \item{\dots}{
additional parameters passed to \code{\link{ddalphaf.train}}
}
}

\value{

  \item{errors}{
  the part of incorrectly classified data (mean)
  }
  \item{errors_sd}{
  the standard deviation of errors
  }
  \item{errors_vec}{
  vector of errors
  }
  \item{time}{
  the mean training time
  }
  \item{time_sd}{
  the standard deviation of training time
  }

}


\seealso{
\code{\link{ddalphaf.train}} to train the functional DD\eqn{\alpha}-classifier, 
\code{\link{ddalphaf.classify}} for classification using functional DD\eqn{\alpha}-classifier, 
\code{\link{ddalphaf.test}} to test the functional DD-classifier on particular learning and testing data,
\code{\link{ddalphaf.getErrorRateCV}} to get error rate of the functional DD-classifier on particular data.
}
\examples{
# load the fdata
df = dataf.growth()

stat <- ddalphaf.getErrorRatePart(dataf = df$dataf, labels = df$labels, 
                          size = 0.3, times = 5,
                          adc.args = list(instance = "avr", 
                                         numFcn = 2, 
                                         numDer = 2))

cat("Classification error rate: ", stat$errors, ".\n", sep = "")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ benchmark }
