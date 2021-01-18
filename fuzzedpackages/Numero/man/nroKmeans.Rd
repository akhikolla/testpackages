\name{nroKmeans}
\alias{nroKmeans}

\title{
K-means clustering
}

\description{
K-means clustering for multi-dimensional data.
}

\usage{
nroKmeans(data, k = 3, subsample = NULL, balance = 0, message = NULL)
}

\arguments{
  \item{data}{
A data frame or a matrix.
  }
  \item{k}{
Number of centroids.
  }  
  \item{subsample}{
Number of randomly selected rows used during a single training cycle.
  }    
  \item{balance}{
Penalty parameter for size difference between clusters.
  }
  \item{message}{
If positive, progress information is printed at the specified
interval in seconds.
  }
}

\details{
The K centroids are determined by Lloyd's algorithm with Euclidean
distances or by using 1 - Pearson correlation as the distance measure.

If \code{subsample} is less than the number of data rows, a random subset of
the specified size is used for each training cycle. By default,
\code{subsample} is set automatically depending on the size of the dataset. 

If \code{balance = 0.0}, the algorithm is applied with no balancing,
if \code{balance = 1.0} all the clusters will be forced to be of equal size.
Intermediate values are permitted. Note that if subsampling is applied,
balancing may become less accurate.
}

\value{
A list with named elements: \code{centroids} is a matrix of the
main results, \code{layout} contains the best-matching centroid labels
and model residuals for each usable data point and \code{history} is the
chronological record of training errors. The subsampling parameter that was
used during training is stored in the element \code{subsample}.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Prepare training data.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- scale.default(dataset[,trvars]) 

# Unbalanced K-means clustering.
km0 <- nroKmeans(data = trdata, k = 5, balance = 0.0)
print(table(km0$layout$BMC))
print(km0$centroids)

# Balanced K-means clustering.
km1 <- nroKmeans(data = trdata, k = 5, balance = 1.0)
print(table(km1$layout$BMC))
print(km1$centroids)
}