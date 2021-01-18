\name{Rcpp_bei-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_bei-class}
\alias{Rcpp_bei}
% \alias{bei}

\title{Class \code{"Rcpp_bei"}}
\description{
This class implements the MCMC sampler for a joint modeling approach to multiple edit-imputation for continuous data. It provides methods for updating and monitoring the sampler. 
}
\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.
}

\keyword{classes}

\section{Fields}{
  \itemize{
    \item \code{Y.input}: input dataset generated from \code{\link{readData}} (replacing NA in \code{Y.original} by -999 and zero values by 0.01).
    \item \code{Y.edited}: current edit-imputed dataset.
    \item \code{K}: number of mixture components (latent classes).
    \item \code{n.occ}: effective number of mixture components.
    % \item \code{Z}: a matrix with the current latent class assignment of each member of the sample.
    \item \code{Prob.A}: ratio of the size of the observed sample to the size of the augmented sample.
    % \item \code{Mu}: Mean vectors of the normal mixture.
    % \item \code{Sigma}: Variance matrix of the normal mixture.
    % \item \code{alpha}: the concentration parameter of the stick breaking prior.
    % \item \code{Pi}: a vector with the mixture weights.
    \item \code{RandomSeed}: random seed.
    \item \code{msg.level}: integer in \{0,1,2\} specifying the level of displayed message; 0: errors only, 1: errors and warnings, 2: all messages. Defaults to 0.
  	\item \code{FaultyRecordID}: record IDs of \code{Y.orig} whose values violate edit rules.
	}
}

\section{Methods}{
  \itemize{
    \item \code{Iterate()}: run a single iteration of MCMC.
    \item \code{Run(iter)}: run \code{iter} iterations of MCMC.
    % \item \code{InitializeSandD()}: function used by \code{\link{createModel}} which finds the initial values of \code{Y.edit}.
    % \item \code{BuildFeasibilityMap()}: function used by \code{\link{createModel}}.
    % \item \code{Simulate_logUnif_case2()}: function used by \code{\link{createModel}}.
  }
}

\details{
Rcpp_bei objects should be created with \code{\link{createModel}}. Please see the example in the demo folder for more detailed explanation. 
}
\references{
Hang J. Kim, Lawrence H. Cox, Alan F. Karr, Jerome P. Reiter and Quanli Wang (2015). "Simultaneous Edit-Imputation for Continuous Microdata", Journal of the American Statistical Association, DOI: 10.1080/01621459.2015.1040881.
}

\examples{
data(SimpleEx)

## read the data
data1 = readData(SimpleEx$D.obs, SimpleEx$Ratio.edit, NULL, 
  SimpleEx$Balance.edit)

## create and initialize the model
# model1 = createModel(data1, K=15)

### run 10 iterations
# model1$Run(10)
# EI_data1 = model1$Y.edited    # store the edit-imputed dataset
}