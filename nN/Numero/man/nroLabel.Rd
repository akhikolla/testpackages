\name{nroLabel}
\alias{nroLabel}

\title{
Label pruning
}

\description{
Optimize the look and selection of labels on map districts.
}

\usage{
nroLabel(topology, values, gap = 2.3)
}

\arguments{
  \item{topology}{
A data frame with K rows and six columns, see details.  
  }
  \item{values}{
A vector of K values or a K x N data frame, where K is the number of map
districts and N is the number of variables.
  }
  \item{gap}{
Minimum distance between map districts with non-empty labels.
  }
}

\details{
The function assigns visible labels for districts based on the absolute
deviations from the average district value. The most extreme districts are
picked first, and then the remaining districts are prioritized based on
their value and distance to the other districts already labeled. Columns
that are listed in the attribute "binary" in \code{values} are given
percentage labels.

Topology can be either the output from \code{\link{nroKohonen}()} or a
data frame in the same format as the element \code{topology} within the
aforementioned output list.
}

\value{
A matrix with K rows and N columns that contains selected labels
for the map districts for each of the columns in \code{values}. The output
has the attribute 'visible' that contains binary flags to guide visibility.
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
sm <- nroTrain(map = sm, data = trdata)

# Assign data points into districts.
matches <- nroMatch(centroids = sm, data = trdata)

# District averages for all variables.
planes <- nroAggregate(topology = sm, districts = matches, data = dataset)

# District labels for cholesterol.
chol <- nroLabel(topology = sm, values = planes[,"CHOL"])
print(head(attr(chol, "visible")))
print(head(chol))

# District labels for all variables.
colrs <- nroLabel(topology = sm, values = planes)
print(head(attr(colrs, "visible")))
print(head(colrs))
}