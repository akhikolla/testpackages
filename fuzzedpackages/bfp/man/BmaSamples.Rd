\name{BmaSamples}
\alias{BmaSamples}

\title{Bayesian model averaging over multiple fractional polynomial
  models} 
\description{
  Draw samples from the Bayesian model average over the models in
  saved in a \code{\link{BayesMfp}}-object. 
}
\usage{
BmaSamples(object, sampleSize = length(object) * 10, postProbs =
posteriors(object), gridList = list(), gridSize = 203, newdata=NULL,
verbose = TRUE, includeZeroSamples=FALSE)  
}

\arguments{
  \item{object}{valid \code{\link{BayesMfp}} object containing the
  models over which to average} 
  \item{sampleSize}{sample size (default is 10 times the number of models)}
  \item{postProbs}{vector of posterior probabilites (will be normalized
  within the function, defaults to the normalized posterior probabilities)} 
  \item{gridList}{optional list of appropriately named grid vectors for FP evaluation,
    default is a length (\code{gridSize} - 2) grid per covariate additional to the observed values
    (two are at the minimum and maximum)}
  \item{gridSize}{see above (default: 203)}
  \item{newdata}{new covariate data.frame with exactly the names (and
    preferably ranges) as before (default: no new covariate data)}  
  \item{verbose}{should information on sampling progress be printed?
    (default)}
  \item{includeZeroSamples}{should the function and coefficient samples
    include zero samples, from models where these covariates are not
    included at all? (default: \code{FALSE}, so the zero samples are
    not included)}
}

\value{
  Return an object of class \code{BmaSamples}, which is a list with
  various elements that describe the \code{BayesMfp} object over which
  was averaged, model frequencies in the samples, the samples themselves
  etc:

  \item{priorSpecs}{the utilized prior specifications}
  \item{termNames}{a list of character vectors containing the names of
    uncertain covariate groups, fractional polynomial terms and fixed
    variables}
  \item{shiftScaleMax}{matrix with 4 columns containing preliminary
    transformation parameters, maximum degrees and cardinalities of the
    powersets of the fractional polynomial terms}
  \item{y}{the response vector}
  \item{x}{the shifted and scaled design matrix for the data}
  \item{randomSeed}{if a seed existed at function call
    (\code{get(".Random.seed", .GlobalEnv)}), it is saved here}
  \item{modelFreqs}{The table of model frequencies in the BMA sample}
  \item{modelData}{data frame containing the normalized posterior
  probabilities of the models in the underlying \code{\link{BayesMfp}}
  object, corresponding log marginal likelihoods, model prior
  probabilities, posterior expected covariance and shrinkage factors,
  coefficients of determination, powers and inclusions, and finally
  model average weights and relative frequencies in the BMA sample.}
\item{sampleSize}{sample size}
\item{sigma2}{BMA samples of the regression variance}
\item{shrinkage}{BMA samples of the shrinkage factor}
\item{fixed}{samples of the intercept}
\item{bfp}{named list of the FP function samples, where each element
  contains one FP covariate and is a matrix (samples x grid), with the
  following attributes:
  \describe{
    \item{whereObsVals}{where in the scaled grid are the originally
      observed covariate values? (integer vector of the indexes)}
    \item{scaledGrid}{numeric vector with the positions of the scaled
      grid points, corresponding to the columns of the samples matrix}
    \item{counter}{how often has this covariate been included in the BMA
      sample? (identical to the number of rows in the samples matrix)}
  }
}
\item{uc}{named list of the uncertain fixed form covariates, where each
  element contains the coefficient samples of one group: in a matrix
  with the attribute \code{counter} as number of samples in the rows,
  and the columns are appropriately named to correspond to the single
  design variables.}
\item{fitted}{fitted values of all models in \code{object}, in a matrix
  with layout models x observations.}
\item{predictions}{samples from the predictive distribution at the
  covariates given in \code{newdata}}
\item{predictMeans}{means of the predictive distribution at the
  covariates given in \code{newdata}}
}

\seealso{\link{BmaSamples Methods}, \code{\link{BayesMfp}}}
\examples{
## construct a BayesMfp object
set.seed(19)

x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)

y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 200, method="exhaustive")

## now draw samples from the Bayesian model average
testBma <- BmaSamples (test)
testBma

## We can also draw predictive samples for new data points, but then
## we need to supply the new data to BmaSamples:
newdata <- data.frame(x1 = rnorm(15),
                      x2 = rbinom(n=15, size=5, prob=0.2) + 1,
                      x3 = rexp(n=15))
testBma <- BmaSamples(test, newdata=newdata)
predict(testBma)

## test that inclusion of zero samples works
testBma <- BmaSamples (test, includeZeroSamples=TRUE)
testBma
}
\keyword{regression}