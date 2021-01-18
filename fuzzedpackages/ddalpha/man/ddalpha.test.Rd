\name{ddalpha.test}
\alias{ddalpha.test}
\title{
Test DD-Classifier
}
\description{
Trains DD-classifier on the learning sequence of the data and tests it on the testing sequence.
}
\usage{
ddalpha.test(learn, test, ...)
}
\arguments{
  \item{learn}{
the learning sequence of the data. Matrix containing training sample where each of \eqn{n} rows is one object of the training sample where first \eqn{d} entries are inputs and the last entry is output (class label).
}
  \item{test}{
the testing sequence. Has the same format as \code{learn}
}
  \item{\dots}{
additional parameters passed to \code{\link{ddalpha.train}}
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
\code{\link{ddalpha.train}} to train the DD-classifier, 
\code{\link{ddalpha.classify}} for classification using DD-classifier, 
\code{\link{ddalpha.getErrorRateCV}} and \code{\link{ddalpha.getErrorRatePart}} to get error rate of the DD-classifier on particular data.
}
\examples{

# Generate a bivariate normal location-shift classification task
# containing 200 training objects and 200 to test with
class1 <- mvrnorm(200, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(200, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
trainIndices <- c(1:100)
testIndices <- c(101:200)
propertyVars <- c(1:2)
classVar <- 3
trainData <- rbind(cbind(class1[trainIndices,], rep(1, 100)), 
                   cbind(class2[trainIndices,], rep(2, 100)))
testData <- rbind(cbind(class1[testIndices,], rep(1, 100)), 
                  cbind(class2[testIndices,], rep(2, 100)))
data <- list(train = trainData, test = testData)

# Train 1st DDalpha-classifier (default settings) 
# and get the classification error rate
stat <- ddalpha.test(data$train, data$test)
cat("1. Classification error rate (defaults): ", 
    stat$error, ".\n", sep = "")

# Train 2nd DDalpha-classifier (zonoid depth, maximum Mahalanobis 
# depth classifier with defaults as outsider treatment) 
# and get the classification error rate
stat2 <- ddalpha.test(data$train, data$test, depth = "zonoid", 
                          outsider.methods = "depth.Mahalanobis")
cat("2. Classification error rate (depth.Mahalanobis): ", 
    stat2$error, ".\n", sep = "")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ benchmark }
