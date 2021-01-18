\name{dknn.train}
\alias{dknn.train}
\title{
Depth-Based kNN
}
\description{
The implementation of the affine-invariant depht-based kNN of Paindaveine and Van Bever (2015). 
}
\usage{
dknn.train(data, kMax = -1, depth = "halfspace", seed = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix containing training sample where each of \eqn{n} rows is one object of the training sample where first \eqn{d} entries are inputs and the last entry is output (class label).
}
  \item{kMax}{
the maximal value for the number of neighbours. If the value is set to -1, the default value is calculated as n/2, but at least 2, at most n-1.
}
  \item{depth}{
Character string determining which depth notion to use; the default value is \code{"halfspace"}. 
Currently the method supports the following depths: "halfspace", "Mahalanobis", "simplicial".
}
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  The returned object contains technical information for classification, including the found optimal value \code{k}.
}
\references{
Paindaveine, D. and Van Bever, G. (2015). Nonparametrically consistent depth-based classifiers. \emph{Bernoulli} \bold{21} 62--82.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dknn.classify}} and \code{\link{dknn.classify.trained}} to classify with the Dknn-classifier.

\code{\link{ddalpha.train}} to train the DD\eqn{\alpha}-classifier.

\code{\link{ddalpha.getErrorRateCV}} and \code{\link{ddalpha.getErrorRatePart}} to get error rate of the Dknn-classifier on particular data (set \code{separator = "Dknn"}).
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

# Train the classifier
# and get the classification error rate
cls <- dknn.train(data$train, kMax = 20, depth = "Mahalanobis")
cls$k
classes1 <- dknn.classify.trained(data$test[,propertyVars], cls)
cat("Classification error rate: ", 
    sum(unlist(classes1) != data$test[,classVar])/200)

# Classify the new data based on the old ones in one step
classes2 <- dknn.classify(data$test[,propertyVars], data$train, k = cls$k, depth = "Mahalanobis")
cat("Classification error rate: ", 
    sum(unlist(classes2) != data$test[,classVar])/200)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ classif }
