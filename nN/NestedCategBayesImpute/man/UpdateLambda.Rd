\name{UpdateLambda}
\alias{UpdateLambda}
\title{
Update lambda.
}
\description{
Update lambda -- the list of matrices of the group-level probabilities for each group-level variable by the group-level latent classes when the weighting/capping option is not used. Each matrix in the list is for each group-level variable.
}
\usage{
UpdateLambda(HHdata_all, G_all, dHH, FF)
}
\arguments{
  \item{HHdata_all}{
  Data matrix for the household-level data from both the original data and the sampled impossible households.
}
  \item{G_all}{
  A vector of the household-level latent class indexes for all households both in the original data and the sampled impossible households.
}
  \item{dHH}{
  A vector containing the number of levels for each household-level variable.
}
  \item{FF}{
  Maximum number of household-level latent classes allowed.
}
}
\details{
Function for obtaining a posterior sample of lambda when the weighting/capping option is not used.
}
\value{
Updated (posterior) value for lambda based on the corresponding full conditional distribution.
}

\author{
Quanli Wang
}

\keyword{ sampler }

