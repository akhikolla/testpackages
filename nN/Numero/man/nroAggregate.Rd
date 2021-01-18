\name{nroAggregate}
\alias{nroAggregate}

\title{
Regional averages on a self-organizing map
}

\description{
Estimate district averages based on assigned map locations for each data point.
}

\usage{
nroAggregate(topology, districts, data = NULL)
}

\arguments{
  \item{topology}{
A data frame with K rows and six columns, see details.  
  }

  \item{districts}{ 
An integer vector of M best-matching districts.
  }

  \item{data}{
A vector of M elements or an M x N matrix of data values. 
  }
}

\details{
Topology can be either the output from \code{\link{nroKohonen}()} or a
data frame in the same format as the element \code{topology} within the
the output from \code{\link{nroKohonen}()}.

The input argument \code{districts} is expected to be the output from
\code{\link{nroMatch}()}.
}

\value{
If the input argument \code{data} is empty, the histogram of the data points
on the map is returned (a vector of K elements).

If data are available, a matrix of K rows and N columns that contains
the average district values after smoothing is returned. The output includes
the attribute 'histogram' that contains data point counts over each data
column. Column names and the attribute 'binary' are copied from the input.
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

# District averages for one variable.
chol <- nroAggregate(topology = sm, districts = matches,
                     data = dataset$CHOL)
print(chol)

# District averages for all variables.
planes <- nroAggregate(topology = sm, districts = matches, data = dataset)
print(head(planes))
}