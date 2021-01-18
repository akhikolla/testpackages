\name{nroKohonen}
\alias{nroKohonen}

\title{
Self-organizing map
}

\description{
Interpolates the initial district profiles of a self-organizing map based
on pre-determined seed profiles.
}

\usage{
nroKohonen(seeds, radius = 3, smoothness = 1.0)
}

\arguments{
  \item{seeds}{
A matrix or a data frame of K rows and N columns.
  }

\item{radius}{
Map radius.
  }

\item{smoothness}{
Rigidity of the map to adapt to regional differences.
  }
}

\value{
A list of named elements: \code{centroids} contains the N-dimensional
district profiles, and \code{topology} is an H x 6 matrix that contains
the 2D spatial layout for the map districts: the first two columns (X, Y)
indicate the positions of districts in Cartesian coordinates, the other
four columns (RADIUS1, RADIUS2, ANGLE1, ANGLE2) define the perimeter of
the district areas for visualisation on a circular map.

Additional parameters are stored as attributes in \code{topology}.

The function is named after Teuvo Kohonen, the inventor of the
self-organizing map.
}

\seealso{
Please see \code{\link{nroKmeans}()} to create the seeds.
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

# Self-organizing map.
sm <- nroKohonen(seeds = km)
print(head(sm$centroids))
print(head(sm$topology))
}