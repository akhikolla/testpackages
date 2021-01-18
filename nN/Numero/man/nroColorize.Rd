\name{nroColorize}
\alias{nroColorize}

\title{
Assign colors based on value
}

\description{
Assign colors to map districts based on the respective district values.
}

\usage{
nroColorize(values, ranges = NULL, amplitudes = 1, palette = "rhodo")
}

\arguments{
  \item{values}{
A vector of K values or a K x N data frame, where K is the number of map
districts and N is the number of variables.
  }
  \item{ranges}{
A data frame with N rows and 2 columns, see details.
  }
  \item{amplitudes}{
Single value or a vector of N elements or a data frame of N rows that
contains the column AMPLITUDE.
  }
  \item{palette}{
One of pre-defined colormap names (see details). 
  }
}

\details{
The argument \code{ranges} sets the minimum and maximum district values
irrespective of the contents of \code{values}. It can be used as a fixed
reference across different colorings to ensure that the same value produces
the same color across function calls.

The argument \code{amplitudes} controls the proportion of the color range
that is available for the district value range. For proportions below 1,
the minimum district value is assigned to a color that is between the
first and middle element in the color palette, and the maximum is assigned
to a color that is between the middle and the last element. If
\code{amplitude} is greater than 1, extreme low and high values are
clipped to the first and last color in the palette, respectively.

Available color palettes include "grey", "fire", "jungle", "miami", "rhodo"
or "tan". Any other word will revert to a rainbow colormap.
}

\value{
A matrix of hexadecimal color codes as strings. The output also includes
the attribute 'contrast' that indicates which colors have a good contrast with
black as opposed to white, and the attribute 'ranges' that contains a copy
of the dynamic ranges across districts.
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

# District colors for cholesterol.
chol <- nroColorize(values = planes[,"CHOL"])
print(head(chol))

# District colors for all variables.
colrs <- nroColorize(values = planes)
print(head(colrs))
}
