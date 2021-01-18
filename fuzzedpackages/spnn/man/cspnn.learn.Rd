\name{cspnn.learn}
\alias{cspnn.learn}
\title{
cspnn.learn
}
\description{
Create or update a Condensed Scale Invariant Probabilistic Neural Network.
}
\usage{
cspnn.learn(set, nn, xr, sigma, category.column = 1)
}
\arguments{
  \item{set}{
data.frame or matrix representing the training set. The first column (default category.column = 1) is used to define the category or class of each observation.
}
  \item{nn}{
(optional) A Condensed Scale Invariant Probabilistic Neural Network object. If provided, the training data set input is concatenated to the current training data set of the neural network. If not provided, a new CSPNN object is created.
}
\item{xr}{
The m by n reference matrix containing optimal parameters for probability estimation. Where m is the number of unique categories and n is the number of input factors used. This matrix must be provided.
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
The function cspnn.learn creates a new Condensed Scale Invariant Probabilistic Neural Network with a given training data set or updates the training data of an existing CSPNN. It sets the parameters: model, set, xr, category.column, categories, sigma, sigmaInverse, k, and n for the CSPNN.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A trained Condensed Scale Invariant Probabilistic Neural Network (CSPNN)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{spnn-package}}, \code{\link{cspnn.predict}}, \code{\link{iris}}
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

# reference matrix must be supplied
# this is not the optimal reference matrix
# this matrix is provided as a simple example
xr <- matrix(c(c(5.00, 3.41, 1.44, 0.24),
               c(5.88, 2.75, 4.23, 1.30),
               c(6.61, 2.97, 5.59, 2.01)),
             nrow = length(unique(trainData$Species)),
             ncol = ncol(trainData) - 1,
             byrow = TRUE)

# fit cspnn
cspnn <- cspnn.learn(set = trainData, xr = xr, category.column = 5)

# estimate probabilities
predictions <- cspnn.predict(nn = cspnn, newData = testData[,1:4])

}
