\name{sampleG}
\alias{sampleG}
\title{
Update household (group) level latent class indexes.
}
\description{
Update household (group) level latent class indexes.
}
\usage{
sampleG(phi, data, omega, pi, ni, HHdata, lambda, Parallel)
}
\arguments{
  \item{phi}{
  Matrix of posterior samples for the individual-level probabilities for each individual-level variable by each pair of group-level and individual-level latent classes.
}
  \item{data}{
  Individual level data.
}
  \item{omega}{
  Matrix of the probabilities for the individual-level latent classes by the group-level latent classes.
}
  \item{pi}{
  Vector of the probabilities for the group-level latent classes.
}
  \item{ni}{
  Vector containing the number of individuals in each household in the data..
}
  \item{HHdata}{
  Household level data.
}
  \item{lambda}{
  A list of matrices of the group-level probabilities for each group-level variable by the group-level latent classes. Each matrix in the list is for each group-level variable.
}
\item{Parallel}{
  Logical indicator for running the function in parallel mode.
}
}
\details{
Function for obtaining a posterior sample of the household-level latent class indexes for all households in the input data based on the corresponding full conditional distribution.
}
\value{
A list with two variables.
  \item{G}{A vector for the updated values of the household-level latent class indexes for all households in the input data.}
  \item{G_Individuals }{The vector \verb{G} expanded to a long format to match the number of individuals in \verb{data}.}
}

\author{
Quanli Wang
}

\keyword{sampler}
\keyword{household level}
