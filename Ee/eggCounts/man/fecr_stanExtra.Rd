\name{fecr_stanExtra}
\alias{fecr_stanExtra}
\title{Model the reduction of faecal egg count using custom models}
\usage{
fecr_stanExtra(preFEC, postFEC, rawCounts = FALSE, preCF = 50, postCF = preCF, 
  modelName = NULL, modelCode = NULL, modelFile = NULL, modelData = NULL,
  nsamples = 2000, nburnin = 1000, thinning = 1, nchain = 2, 
  ncore = 1, adaptDelta = 0.95, verbose = FALSE)
}
\arguments{
\item{preFEC}{numeric vector. Pre-treatment faecal egg counts. Not required if \code{modelData} is supplied.}

\item{postFEC}{numeric vector. Post-treatment faecal egg counts. Not required if \code{modelData} is supplied.}

\item{rawCounts}{logical. If TRUE, \code{preFEC} and \code{postFEC} correspond to raw counts
(as counted on equipment). Otherwise they correspond to calculated epgs (raw counts times correction factor).
Defaults to \code{FALSE}. Not required if \code{modelCode} or \code{modelFile} is supplied.}

\item{preCF}{a positive integer or a vector of positive integers. Pre-treatment correction factor(s). Not required if \code{modelCode} or \code{modelFile} is supplied.}

\item{postCF}{a positive integer or a vector of positive integers. Post-treatment correction factor(s). Not required if \code{modelCode} or \code{modelFile} is supplied.}

\item{modelName}{string. One of four availale models ("Po", "UPo", "ZIPo", "ZIUPo") from \strong{eggCountsExtra} package, which corresponds to outlier-adjusted version of paired, unpaired, paired with zero inflation and unpaired with zero inflation models. Not required if \code{modelCode} or \code{modelFile} is supplied.}

\item{modelCode}{stan model code. Not required when \code{modelName} or \code{modelFile} is supplied.}

\item{modelFile}{stan model file with file extension \sQuote{*.stan}. Not required when \code{modelName} or \code{modelCode} is supplied.}

\item{modelData}{stan data list. A named list or environment providing the data for the model, or a character vector for all the names of objects in the current enviroment used as data. Not required when \code{modelName} is supplied.}

\item{nsamples}{a positive integer. Number of samples for each chain (including burn-in samples).}

\item{nburnin}{a positive integer. Number of burn-in samples.}

\item{thinning}{a positive integer. Thinning parameter, i.e. the period for saving samples.}

\item{nchain}{a positive integer. Number of chains.}

\item{ncore}{a positive integer. Number of cores to use when executing the chains in parallel.}

\item{adaptDelta}{numeric. The target acceptance rate, a numeric value between 0 and 1.}

\item{verbose}{logical. If TRUE, prints progress and debugging information.}

}
\value{
Prints out the posterior summary of \code{FECR} as the reduction, \code{meanEPG.untreated} as the mean pre-treatment epg, and \code{meanEPG.treated} as the mean after-treatment epg. The posterior summary contains the mean, standard deviation (sd), 2.5\%,  50\% and 97.5\% percentiles, the 95\% highest posterior density interval (HPDLow95 and HPDHigh95) and the posterior mode. 

The returned value is a list that consists of:
\item{stan.model}{an object of class \code{\link[rstan]{stanmodel-class}} that was used}
\item{stan.samples}{an object of S4 class \code{\link[rstan]{stanfit}} representing the fitted results}
\item{posterior.summary}{a data.frame that is the same as the printed posterior summary. Not available for custom models.}
}
\description{
Models the reduction in faecal egg counts with custom model formulation using Stan modelling language (for advanced users).
}
\details{
If  \code{modelName} is one of c("Po", "UPo", "ZIPo", "ZIUPo"), then outlier-adjusted models are used. 

\itemize{
\item In paired models, outliers are those counts with postFEC > preFEC. Outlier weights are assigned as the inverse of postFEC/preFEC,
\item In unpaired models, outliers are those counts with postFEC greater than the 95th percentile of a Poisson distribution, where the Poisson mean is computed based on the mean of postFEC excluding postFEC > Q3 + 1.5*IQR. Q3 is the 75th percentile and IQR is the interquartile range. The lowest outlier weight is assigned as 0.01, and other outliers assigned proportionally.
\item In both cases, non-outliers are assigned with outlier weight = 1.
} 

The first time each model is applied, it can take up to 20 seconds for Stan to compile the model. 

The default number of samples per chain is 2000, with 1000 burn-in samples. Normally this is sufficient in Stan. If the chains do not converge, one should tune the MCMC parameters until convergence is reached to ensure reliable results.
}

\author{
Craig Wang
}

\examples{
\dontrun{
library(eggCountsExtra)
data(epgs) ## load sample data

## apply paired model with outliers 
model1 <- fecr_stanExtra(epgs$before, epgs$after, rawCounts=FALSE, 
         preCF=10, modelName = "Po")
samples <- stan2mcmc(model1$stan.samples)
fecr_probs(model1$stan.samples, threshold = 0.99)

## apply a simple custom model
code <- "data{
  int J; // number of animals
  int y_before[J]; // after treatment McMaster count
  int y_after[J]; // before treatment McMaster count
}
parameters{
  real<lower=0> mu;
  real<lower=0,upper=1> delta;
}
model{
  mu ~ gamma(1,0.001);
  delta ~ beta(1,1);
  y_before ~ poisson(mu);
  y_after ~ poisson(mu*delta);
}"

dat <- list(J = nrow(epgs), y_before = epgs$before,
            y_after = epgs$after)
model2 <- fecr_stanExtra(modelCode = code, modelData = dat)
}}
\keyword{models}

