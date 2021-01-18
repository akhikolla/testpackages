\name{nroPermute}
\alias{nroPermute}

\title{
Permutation analysis of map layout
}

\description{
Estimate the dynamic range and statistical significance for
regional patterns on a self-organizing maps using permutations.
}

\usage{
nroPermute(map, districts, data, n = 1000, message = NULL)
}

\arguments{
  \item{map}{
A list object in the format from \code{\link{nroTrain}()}.
  }
  
  \item{districts}{
An integer vector of M best matching districts.
  }
  
  \item{data}{
A numeric vector of M values or an M x N matrix (or data frame), where
M is the number of data points and N is the number of variables.
  }
  
  \item{n}{
Maximum number of permutations per variable.
  }

  \item{message}{
If positive, progress information is printed at the specified
interval in seconds.
  }
}

\details{
The input argument \code{map} must contain the map topology and the
centroid profiles as returned by the functions \code{\link{nroKmeans}()},
\code{\link{nroKohonen}()}, or \code{\link{nroTrain}()}.

The input argument \code{districts} must contain integers between 1 and K,
where K is the number map units. Any other values will be ignored.

Training variables and data points are detected by the column names of
\code{map$centroids}, the attribute "variables" in \code{districts} and
the names of elements in \code{districts}.
}

\value{
A data frame with eight columns: P.z is a parametric estimate for statistical
significance, P.freq is the frequency-based estimate for statistical
signicance, and Z is the estimated z-score of how far the
observed map plane is from the average randomly generated layout.
N.data indicates how many data values were used and N.cycles tells the
number of completed permutations. AMPLITUDE is a dynamic range modifier
for colors that can be used in \code{\link{nroColorize}()}.

The output also contains the attribute 'zbase' that indicates the
normalization factor for the color amplitudes.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Set row names.
rownames(dataset) <- paste("r", 1:nrow(dataset), sep="")

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

# Estimate statistics for cholesterol
chol <- nroPermute(map = sm, districts = matches, data = dataset$CHOL)
print(chol[,c("TRAINING", "Z", "P.z", "P.freq")])

# Estimate statistics.
stats <- nroPermute(map = sm, districts = matches, data = dataset)
print(stats[,c("TRAINING", "Z", "P.z", "P.freq")])
}