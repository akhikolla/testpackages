\name{getStuffForDIC}
\alias{getStuffForDIC}
\title{
Compute information criteria
}
\description{
This function computes four information criteria for a given run of the \code{fabMix} algorithm, namely: AIC, BIC, DIC and \eqn{\mbox{DIC}_2}. Given various runs with different number of factors, the selected model corresponds to the one with the smalled value of the selected criterion.
}
\usage{
getStuffForDIC(sameSigma, sameLambda, isotropic, x_data, outputFolder, q, burn, 
				Km, normalize, discardLower)
}
\arguments{
  \item{sameSigma}{
	Logical value indicating whether the parameterization with the same variance of errors per component is used. Default: TRUE.
}
  \item{sameLambda}{
	Logical value indicating whether the parameterization with same loadings per component is used. Default: FALSE.
}
  \item{isotropic}{
	Logical value indicating whether the parameterization with isotropic error variance per component is used. Default: FALSE.
}
  \item{x_data}{
	Observed data.
}
  \item{outputFolder}{
	Name of the folder where the \code{fabMix} function has saved its output.
}
  \item{q}{
	Number of factors. Note that this should coincide with the number of factors in the \code{fabMix} run.
}
  \item{burn}{
	Discard observations as burn-in period (optional).
}
  \item{Km}{
	Number of components in the overfitted mixture model. Note that this should coincide with the same entry in the \code{fabMix} run.
}
  \item{normalize}{
	Should the observed data be normalized? Note that this should coincide with the same entry in the \code{fabMix} run. Default value: TRUE.
}
  \item{discardLower}{
	Discard draws with log-likelihood values lower than the specific quantile. This applied only for the DIC computation. 
}
}
\details{
If necessary, more details than the description above
}
\value{
	The information criteria are saved to the \code{informationCriteria_map_model.txt} file in the code{outputFolder}. 
}
\author{
Panagiotis Papastamoulis
}
\note{
It is well known that DIC tends to overfit, so it advised to compare models with different number of factors using AIC or BIC. The main function of the package uses BIC. 
}

