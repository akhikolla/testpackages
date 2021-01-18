\name{nroMatch}
\alias{nroMatch}

\title{
Best-matching districts
}

\description{
Compare multi-dimensional data points against the district
profiles of a self-organizing map (SOM).
}

\usage{
nroMatch(centroids, data)
}

\arguments{
  \item{centroids}{
Either a matrix, a data frame or a list that contains the element
\code{centroids}.
  }
  
  \item{data}{
A data matrix with identical column names to the centroid matrix.
  }
}

\details{
The input argument \code{centroids} can be a matrix or a data frame that
contains multivariable data profiles organized row-wise. It can also be
the output list object from \code{\link{nroKmeans}()} or
\code{\link{nroTrain}()}.
}

\value{
A vector of integers with elements corresponding to the rows in
\code{data}. Each element contains the index of the best matching
row from \code{centroids}.

The vector also has the attribute 'quality' that contains three columns:
RESIDUAL is the distance between a point and a centroid in data space
(shorter is better), RESIDUAL.z is a scale-independent version of RESIDUAL
if the mean residual and standard deviation are available from training
history, and COVERAGE shows the proportion of data elements that were
available for matching.

The names of the columns that were used for matching are stored in the
attribute \code{variables}.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Prepare training data.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- scale.default(dataset[,trvars]) 

# K-means clustering.
km <- nroKmeans(data = trdata, k = 10)

# Assign data points into districts.
matches <- nroMatch(centroids = km, data = trdata)
print(head(attr(matches,"quality")))
print(table(matches))
}