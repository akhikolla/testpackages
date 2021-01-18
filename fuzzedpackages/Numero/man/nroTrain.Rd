\name{nroTrain}
\alias{nroTrain}
\title{
Train self-organizing map
}
\description{
Iterative algorithm to adapt a self-organizing map (SOM) to a set of
multivariable data.
}
\usage{
nroTrain(map, data, subsample = NULL, message = NULL)
}
\arguments{
  \item{map}{
A list object as returned by \code{\link{nroKohonen}()}.
  }
  \item{data}{
A matrix or a data frame.
  }
  \item{subsample}{
Number of rows used during a single training cycle.
  }
  \item{message}{
If positive, progress information is printed at the specified
interval in seconds.
  }
}

\details{
The map is fitted according to columns that are found both in the
SOM centroids and the input data.

If \code{subsample} is less than the number of data rows, a random subset of
the specified size is used for each training cycle. By default,
\code{subsample} is set automatically depending on the size of the dataset. 
}

\value{
A copy of the list object \code{map}, where the element \code{centroids} is
updated according to the data patterns. The quantization errors during
training are stored in the element \code{history}. The subsampling
parameter that was used during training is stored in the element
\code{subsample}.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Prepare training data.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- scale.default(dataset[,trvars]) 

# K-means clustering.
km <- nroKmeans(data = trdata)

# Train with full data.
sm <- nroKohonen(seeds = km)
sm <- nroTrain(map = sm, data = trdata, subsample = nrow(trdata))
print(sm$history)

# Train with subsampling.
sm <- nroKohonen(seeds = km)
sm <- nroTrain(map = sm, data = trdata, subsample = 200)
print(sm$history)
}
