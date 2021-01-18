\name{spnn.learn}
\alias{spnn.learn}
\title{
spnn.learn
}
\description{
Create or update a Scale Invariant Probabilistic Neural Network.
}
\usage{
spnn.learn(set, nn, sigma, category.column = 1)
}
\arguments{
  \item{set}{
data.frame or matrix representing the training set. The first column (default category.column = 1) is used to define the category or class of each observation.
}
  \item{nn}{
(optional) A Scale Invariant Probabilistic Neural Network object. If provided, the training data set input is concatenated to the current training data set of the neural network. If not provided, a new SPNN object is created.
}
  \item{sigma}{
An n by n square matrix of smoothing parameters where n is the number of input factors. Defaults to using the covariance matrix of the training data set excluding the category.column.
}
  \item{category.column}{
The column number of category data. Default is 1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function spnn.learn creates a new Scale Invariant Probabilistic Neural Network with a given training data set or updates the training data of an existing SPNN. It sets the parameters: model, set, category.column, categories, sigma, sigmaInverse, k, and n for the SPNN.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A trained Scale Invariant Probabilistic Neural Network (SPNN)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{spnn-package}}, \code{\link{spnn.predict}}, \code{\link{iris}}
}
\examples{
library(spnn)
library(datasets)

data(iris)

# shuffle the iris data set
indexRandom <- sample(1:nrow(iris), size = nrow(iris), replace = FALSE)

# use 100 observations for training set
trainData <- iris[indexRandom[1:100],]

# use remaining observations for testing
testData <- iris[indexRandom[101:length(indexRandom)],]

# fit spnn
spnn <- spnn.learn(set = trainData, category.column = 5)

# estimate probabilities
predictions <- spnn.predict(nn = spnn, newData = testData[,1:4])

}
