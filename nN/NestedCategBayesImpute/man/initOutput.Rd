\name{initOutput}
\alias{initOutput}

\title{
Set the output structure for saving posterior samples of parameters.
}
\description{
Set the output structure for saving posterior samples of parameters.
}
\usage{
initOutput(data, hyper, mc)
}

\arguments{
  \item{data}{
  A list object including all the necessary data variables needed by the sampler.; output of the \verb{initData} function.
}
  \item{hyper}{
  Hyper parameters for priors.
}
  \item{mc}{
  MCMC parameters.
}
}
%\details{
% Currently, the program is harded code to the specific data format used in the the paper. Needs to be generalized once a standard foramt is %defined and the general format to define structural zeros is finalized.
%}
\value{
A list of output parameters to be saved.
\item{alphaout}{Vector of posterior samples for the concentration parameter in the Dirichlet process for the group-level latent classes.}
\item{betaout}{Vector of posterior samples for the concentration parameter in the Dirichlet process for the individual-level latent classes. Currently, this is assumed to be the same within all group-level classes.}
\item{piout}{Matrix of posterior samples for the vector of probabilities for the group-level latent classes.}
\item{omegaout}{3D array of posterior samples for the matrix of probabilities for the group-level and individual-level latent class pairs.}
\item{nout}{Vector of posterior samples for the total number of impossible households sampled.}
\item{extrasize}{Matrix of posterior samples for the number of impossible households sampled, split by household size.}
\item{F_occupied}{Vector of posterior samples for the number of occupied household-level latent classes.}
\item{S_occupied_max}{Vector of posterior samples for the max number of occupied individual-level latent classes.}
\item{elapsed_time}{Vector of time taken to run each iteration.}
\item{newphiout}{3D array of posterior samples for the individual-level probabilities for each individual-level variable by each pair of group-level and individual-level latent classes.}
\item{lambdaout}{A list of an array of posterior samples for the group-level probabilities for each group-level variable. Each array in the list is for each group-level variable.}
}

\author{
Quanli Wang, Olanrewaju Akande
}

