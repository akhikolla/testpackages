\name{initMissing}
\alias{initMissing}
\title{
Initilize the misising data structure from input data
}
\description{
Initilize the misising data structure from input data
}
\usage{
initMissing(data,struc_zero_variables,miss_batch)
}
\arguments{
  \item{data}{A list that holds all input data info.}
  \item{struc_zero_variables}{column indexes for the variables that define structural zeros like age and relate (including those for the household head).}
  \item{miss_batch}{initial number of batches to sample for each household with missing data.}
}
