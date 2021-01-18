\name{nroDestratify}
\alias{nroDestratify}

\title{Mitigate data stratification}

\description{Removes differences in value distribution within
subsets of data points.}

\usage{
nroDestratify(data, labels)
}

\arguments{
  \item{data}{
A matrix or a data frame with M rows.
  }

  \item{labels}{
A vector of M subset labels.
  }
}

\details{
Only non-binary numerical columns are processed, the rest will be
excluded from the results.

The de-stratification algorithm is based on ranked data: the distribution
of each subset will be mapped to the pooled distribution over all subsets
by matching subset-specific ranking with ranking of all values.
}

\value{
A matrix of de-stratified values. The output also includes the attribute
'incomplete' that lists those columns where (some of) the values were set
to missing due to processing failures.
}

\examples{
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname)

# Remove sex differences for creatinine.
creat <- nroDestratify(dataset$CREAT, dataset$MALE)

# Compare creatinine distributions.
men <- which(dataset$MALE == 1)
women <- which(dataset$MALE == 0)
print(summary(dataset[men,"CREAT"]))
print(summary(dataset[women,"CREAT"]))
print(summary(creat[men]))
print(summary(creat[women]))

# Remove sex differences (produces warnings for binary traits).
ds <- nroDestratify(dataset, dataset$MALE)

# Compare HDL2C distributions.
print(summary(dataset[men,"HDL2C"]))
print(summary(dataset[women,"HDL2C"]))
print(summary(ds[men,"HDL2C"]))
print(summary(ds[women,"HDL2C"]))
}