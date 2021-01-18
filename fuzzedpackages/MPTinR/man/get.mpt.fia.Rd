\name{get.mpt.fia}
\alias{get.mpt.fia}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenient function to get FIA for MPT
}
\description{
\code{get.mpt.fia} is a comfortable wrapper for the R-port of Wu, Myung, and Batchelder's (2010) BMPTFIA \code{\link{bmpt.fia}}. It takes data, a model file, and (optionally) a restrictions file, computes the context-free representation of this file and then feeds this into \code{\link{bmpt.fia}} which returns the FIA.
}
\usage{
get.mpt.fia(data, model.filename, restrictions.filename = NULL, Sample = 2e+05, 
	model.type = c("easy", "eqn", "eqn2"), round.digit = 6, 
	multicore = FALSE, split = NULL, mConst = NULL)
}
\arguments{
  \item{data}{
Same as in \code{\link{fit.mpt}}
  }
  \item{model.filename}{
Same as in \code{\link{fit.mpt}}
}
  \item{restrictions.filename}{
Same as in \code{\link{fit.mpt}}
}
  \item{Sample}{
The number of random samples to be drawn in the Monte Carlo algorithm. Default is 200000.
}
  \item{model.type}{
Same as in \code{\link{fit.mpt}}
}
  \item{round.digit}{
scalar numeric indicating to which decimal the ratios between ns in trees should be rounded (for minimizing computations with differing ns, see Details)
}
  \item{multicore}{
Same as in \code{\link{bmpt.fia}}
}
  \item{split}{
Same as in \code{\link{bmpt.fia}}
}
  \item{mConst}{
A constant which is added in the Monte Carlo integration to avoid numerical underflows and is later subtracted (after appropriate transformation). Should be a power of 2 to avoid unnecessary numerical imprecision. 
}
}
\details{
This function is called from \code{\link{fit.mpt}} to obtain the FIA, but can also be called independently.

It performs the following steps:\cr
1.) Equality restrictions (if present) are applied to the model.\cr
2.) The representation of the model as equations is transformed to the string representation of the model into the context-free language of MPT models (L-BMPT; Purdy & Batchelder, 2009). For this step to be successful it is \emph{absolutely necessary} that the equations representing the model perfectly map the tree structure of the MPT. That is, the model file is only allowed to contain parameters, their negations (e.g., \code{Dn} and \code{(1 - Dn)}) and the operators + and *, but nothing else. Simplifications of the equations will seriously distort this step. This step is achieved by calling \code{\link{make.mpt.cf}}.\cr
Note that inequality restrictions are not included in this transformation.\cr
3.) The context free representation of the model is then fed into the MCMC function computing the FIA (the port of BMPTFIA provided by Wu, Myung, & Batchelder, 2010; see \code{\link{bmpt.fia}}).\cr
If inequality restrictions are present, these are specified in the call to \code{\link{bmpt.fia}}.

For multi-individual data sets (i.e., \code{data} is a \code{matrix} or \code{data.frame}), \code{get.mpt.fia} tries to minimize computation time. That is done by comparing the ratios of the number items between trees. To not run into problems related to floating point precision, these values are rounded to \code{round.digit}. Then, \code{get.mpt.fia} will only call \code{\link{bmpt.fia}} as many times as there are differing ratios. As a consequence, the final penalty factor for FIA (\code{CFIA}) is calculated by \code{get.mpt.fia}, \emph{without providing confidence intervals for the penalty factor}.
}
\value{
A \code{data.frame} containing the results as returned by \code{\link{bmpt.fia}}:
\item{CFIA}{The FIA complexity value of the model with the corresponding confidence interval \code{CI.l} (lower bound) and \code{CI.u} (upper bound).}
\item{lnInt}{The log integral term in C_FIA (Wu, Myung, & Batchelder, 2010a; Equation 7) for models without inequality constraints. When inequality constraints are present, 'lnInt' does not take into account the change in the normalizing constant in the proposal distribution and must be adjusted with the output argument `lnconst'. The corresponding confidence interval ranges from \code{CI.lnint.l} (lower bound) to \code{CI.lnint.u} (upper bound).}
\item{lnconst}{When inequality constraints are present, \code{lnconst} serves as an adjustment of code{lnInt}. It estimates the logarithm of the proportion of parameter space [0,1]^S that satisfies those inequality constraints, and the log integral term is given by \code{lnInt+lnconst}.

The next [two] output argument [CI.lnconst] give the Monte Carlo confidence interval of `lnconst'. [.l = lower & .u = upper bound of the CI]}
}
\references{
Purdy, B. P., & Batchelder, W. H. (2009). A context-free language for binary multinomial processing tree models. \emph{Journal of Mathematical Psychology}, 53, 547-561.

Wu, H., Myung, J.I., & Batchelder, W.H. (2010). Minimum description length model selection of multinomial processing tree models. \emph{Psychonomic Bulletin & Review}, 17, 275-286.
}
\author{
Henrik Singmann
}

\seealso{
calls \code{\link{bmpt.fia}}\cr
is called by \code{\link{fit.mpt}}, the main function for fitting MPT models
}
\examples{

# Get the FIA for the 40 datasets from Broeder & Schuetz (2009, Experiment 3)
# for the 2HTM model with inequality restrictions
# (Can take a while.)

data(d.broeder)
m.2htm <- system.file("extdata", "5points.2htm.model", package = "MPTinR")
i.2htm <- system.file("extdata", "broeder.2htm.ineq", package = "MPTinR")

get.mpt.fia(d.broeder, m.2htm, Sample = 1000) # Way too little samples
get.mpt.fia(d.broeder, m.2htm, i.2htm, Sample = 1000)

\dontrun{
# should produce very similar results:
get.mpt.fia(d.broeder, m.2htm, i.2htm)
get.mpt.fia(d.broeder, m.2htm, i.2htm, mConst = 2L^8)
}
  }



