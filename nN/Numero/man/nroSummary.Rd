\name{nroSummary}
\alias{nroSummary}

\title{
Estimate subgroup statistics
}

\description{
Combine data points that reside in districts that belong to a larger region
into a subgroup; compare descriptive statistics between subgroups.
}

\usage{
nroSummary(data, districts, regions = NULL, categlim = 8, capacity = 10)
}

\arguments{
  \item{data}{
A vector of named M elements or an M x N matrix of data values with row names. 
  }
  
  \item{districts}{ 
An integer vector of M named elements that indicate the best match out of
K districts for each row name in the data matrix, please see \link{nroMatch}
for an example.
  }

  \item{regions}{
An vector of K elements or a data frame of K rows that defines if a
district belongs to a larger region (i.e. a subgroup), see details.
  }
  
  \item{categlim}{
The threshold for the number of unique values before a variable is
considered continuous.
  }

  \item{capacity}{
Maximum number of subgroups to compare.
  }
}

\details{
If defined, the region vector should have K elements where K is the total
number of map districts.

The region input can also be a data frame of K rows where the column
REGION will be used for assigning district to regions, and REGION.label
will be used as the character label as seen on the map, see the output
from \code{\link{nroPlot}()} for an example.

Districts and data points are connected by comparing element names in
\code{districts} and names or row names of \code{data}.

Districts and regions are connected by comparing element values in
\code{districts} and names or row names of \code{regions}.

If the region vector is empty, each district is automatically assigned to
its own region.

Safeguards are in place to prevent crashes from empty categories; this
reduces statistical power slightly when numbers are small.
}

\value{
A data frame of summary statistics that contains a row for every
combination of subgroups and variables. The chi-squared test is used for
comparisons with respect to categorical variables, and rank-regulated t-test
and ANOVA are applied to continuous variables. Region labels for each row
are stored in the attribute 'labels' and a list that contains the subsets
of rows in each region is stored in the attribute 'subgroups'.
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

# Calculate district averages for urinary albumin.
plane <- nroAggregate(topology = sm, districts = matches,
                      data = dataset$uALB)
plane <- as.vector(plane)

# Assign subgroups based on urinary albumin.
regns <- rep("HighAlb", length.out=length(plane))
regns[which(plane < quantile(plane, 0.67))] <- "MiddleAlb"
regns[which(plane < quantile(plane, 0.33))] <- "LowAlb"

# Add label info and make a data frame.
regns <- data.frame(REGION=regns, REGION.label="",
    stringsAsFactors=FALSE)
regns[which(regns$REGION == "HighAlb"),"REGION.label"] <- "H"
regns[which(regns$REGION == "MiddleAlb"),"REGION.label"] <- "M"
regns[which(regns$REGION == "LowAlb"),"REGION.label"] <- "L"

# Calculate summary statistics.
st <- nroSummary(data = dataset, districts = matches, regions = regns)
print(st[,c("VARIABLE","SUBGROUP","MEAN","P.chisq","P.t","P.anova")])
}