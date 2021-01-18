\name{numero.create}
\alias{numero.create}

\title{Create a self-organizing map}

\description{Set up a self-organizing map and train it with data}

\usage{
numero.create(data, radius = NULL, smoothness = NULL, subsample = NULL)
}

\arguments{
  \item{data}{
A matrix or a data frame. 
  }
  \item{radius}{
Map radius.
  }
\item{smoothness}{
Rigidity of the map to adapt to regional differences.
  }
  \item{subsample}{
Number of data points used during a single training cycle.
  }
}

\details{
The parameter \code{subsample} sets the number of data points that are
randomly picked for each training cycle; if the number is substantially less
than the size of the dataset, the function will finish quicker.
}

\value{
A list with named elements: \code{data} contains the training data,
\code{kmeans} is the output from \code{\link{nroKmeans}()} during the
initialiation of the SOM, \code{map} is the finished self-organising map
from \code{\link{nroTrain}()} and \code{layout} contains the output from
\code{\link{nroMatch}()} for the training data points.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Set identities and manage missing data.
dataset <- numero.clean(dataset, identity = "INDEX")

# Prepare training set.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
trdata <- numero.prepare(data = dataset, variables = trvars)

# Create a self-organizing map.
modl <- numero.create(data = trdata)
}