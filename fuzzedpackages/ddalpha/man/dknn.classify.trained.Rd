\name{dknn.classify.trained}
\alias{dknn.classify.trained}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Depth-Based kNN
}
\description{
The implementation of the affine-invariant depth-based kNN of Paindaveine and Van Bever (2015). 
}
\usage{
dknn.classify.trained(objects, dknn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{objects}{
Matrix containing objects to be classified; each row is one \eqn{d}-dimensional object.
}
  \item{dknn}{
Dknn-classifier (obtained by \code{\link{dknn.train}}).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
List containing class labels, or character string "Ignored" for the outsiders if "Ignore" was specified as the outsider treating method.
}
\references{
Paindaveine, D. and Van Bever, G. (2015). Nonparametrically consistent depth-based classifiers. \emph{Bernoulli} \bold{21} 62--82.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dknn.train}} to train the Dknn-classifier.

\code{\link{dknn.classify}} to classify with the Dknn-classifier.

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