\name{dealWithLabelSwitching}
\alias{dealWithLabelSwitching}
\title{
	Apply label switching algorithms
}
\description{
	This functions is a wrapper for the \code{label.switching} package and applies the \code{ECR} and \code{ECR.ITERATIVE.1} algorithms. The model may have the same variance of error terms per cluster or not. 
}
\usage{
dealWithLabelSwitching(sameSigma, x_data, outputFolder, q, burn, 
		z.true, compute_regularized_expression, Km)
}
\arguments{
  \item{sameSigma}{
	Logical value indicating whether the parameterization with the same error variance per cluster is used.
}
  \item{x_data}{
	Data
}
  \item{outputFolder}{
	Name of the folder where the \code{fabMix} function has saved its output
}
  \item{q}{
	Number of factors
}
  \item{burn}{
	Discard observations as burn-in period (optional). 
}
  \item{z.true}{
	An (optional) vector of cluster assignments which is considered as the groun-truth clustering of the data. Useful for direct comparisons against the real parameter values in simulated data. 
}
  \item{compute_regularized_expression}{
	Logical. Should regularized expression be computed?
}
  \item{Km}{
	Number of components in the overfitted mixture model.
}
}
\value{
	The following files are produced in the output folder:
}
\references{
Papastamoulis, P. (2016). \code{label.switching}: An R Package for Dealing with the Label Switching Problem in MCMC Outputs. Journal of Statistical Software, 69(1), 1-24.
}
\author{
Panagiotis Papastamoulis
}


