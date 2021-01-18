\name{ddalphaf.getErrorRateCV}
\alias{ddalphaf.getErrorRateCV}
\title{
Test Functional DD-Classifier
}
\description{
Performs a cross-validation procedure over the given data. 
On each step every \code{numchunks} observation is removed from the data, the functional DD-classifier is trained on these data and tested on the removed observations.
}
\usage{
ddalphaf.getErrorRateCV (dataf, labels, numchunks = 10, disc.type = c("LS", "comp"),  ...)
}
\arguments{
  \item{dataf}{
list containing lists (functions) of two vectors of equal length, named "args" and "vals": arguments sorted in ascending order and corresponding them values respectively
}
  \item{labels}{
list of output labels of the functional observations
}
  \item{numchunks}{
number of subsets of testing data. Equals to the number of times the classifier is trained.
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
  the part of incorrectly classified data
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
\code{\link{ddalphaf.getErrorRatePart}} to perform a benchmark study of the functional DD-classifier on particular data.
}
\examples{
# load the fdata
df = dataf.growth()

stat <- ddalphaf.getErrorRateCV(dataf = df$dataf, labels = df$labels, 
                                numchunks = 5,
                                adc.args = list(instance = "avr", 
                                                numFcn = 2, 
                                                numDer = 2))
                                                
cat("Classification error rate: ", stat$errors, ".\n", sep = "")

    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ benchmark }
