\name{BayesMfp}
\alias{BayesMfp}
\alias{bfp}
\alias{uc}

\title{Bayesian model inference for multiple fractional polynomial models}
\description{
  Bayesian model inference for multiple fractional polynomial
  models is conducted by means of either exhaustive model space
  evaluation or posterior model sampling.
}
\usage{
BayesMfp(formula = formula(data), data = parent.frame(), family =
gaussian, priorSpecs = list(a = 4, modelPrior = "flat"), method =
c("ask", "exhaustive", "sampling"), subset = NULL, na.action = na.omit,
verbose = TRUE, nModels = NULL, nCache=1e9L, chainlength = 1e5L)

bfp(x, max = 2, scale = TRUE, rangeVals=NULL)

uc(x)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{model formula}
  \item{data}{optional data.frame for model variables (defaults to the
    parent frame)}
  \item{family}{distribution and link: only gaussian("identity") supported at the moment}
  \item{priorSpecs}{prior specifications, see details}
  \item{method}{which method should be used to explore the  posterior
    model space? (default: ask the user)}
  \item{subset}{optional subset expression}
  \item{na.action}{default is to skip rows with missing data, and no other
    option supported at the moment}
  \item{verbose}{should information on computation progress be given? (default)}
  \item{nModels}{how many best models should be saved? (default: 1\% of
    the explored models or the chainlength, 1 would mean only the
    maximum a posteriori [MAP] model)}
  \item{nCache}{maximum number of best models to be cached at the same
    time during the model sampling (only has an effect if sampling has
    been chosen as method)}  
  \item{chainlength}{length of the model sampling chain (only has an
    effect if sampling has been chosen as method)}
  \item{x}{variable}
  \item{max}{maximum degree for this FP (default: 2)}
  \item{scale}{use pre-transformation scaling to avoid numerical
    problems? (default)}
  \item{rangeVals}{extra numbers if the scaling should consider values
    in this range. Use this argument if you have test data with larger
    range than the training range.}   
  }
\details{
  The formula is of the form
  \code{y ~ bfp (x1, max = 4) + uc (x2 + x3)}, that is, the
  auxiliary functions \code{\link{bfp}} and \code{\link{uc}} must be
  used for defining the fractional polynomial and uncertain fixed form
  covariates terms, respectively. There must be an intercept, and no
  other fixed covariates are allowed. All \code{max} arguments of the
  \code{\link{bfp}} terms must be identical.

  The prior specifications are a list:
  \describe{
    \item{a}{hyperparameter for hyper-g prior which must be greater than
      3 and is recommended to be not greater than 4 (default is 4)}
    \item{modelPrior}{choose if a flat model prior (default,
      \code{"flat"}), a model prior favoring 
      sparse models explicitly (\code{"sparse"}), or a dependent model
    prior (\code{"dependent"}) should be used.}
  }

  If \code{method = "ask"}, the user is prompted with the maximum
  cardinality of the model space and can then decide whether to use
  posterior sampling or the exhaustive model space evaluation.

  Note that if you specify only one FP term, the exhaustive model search
  must be done, due to the structure of the model sampling algorithm.
  However, in reality this will not be a problem as the model space will
  typically be very small.
}
\value{
  Returns an object of class \code{BayesMfp} that inherits from list. It
  is essentially a list of models. Each model is a list and has the
  following components: 

  \item{powers}{a list of numeric vectors, where each vector contains
    the powers of the covariate that its name denotes.}
  \item{ucTerms}{an integer vector of the indices of uncertain fixed
    form covariates that are present in the model.}
  \item{logM}{log marginal likelihood}
  \item{logP}{log prior probability}
  \item{posterior}{normalized posterior probability, and if model
    sampling was done, the frequency of the model in the sampling
    algorithm} 
  \item{postExpectedg}{posterior expected covariance factor g}
  \item{postExpectedShrinkage}{posterior expected shrinkage factor
    t=g/(g + 1)}
  \item{R2}{usual coefficient of determination for the linear model}

  Subsetting the object
  with \code{\link{[.BayesMfp}} returns again a \code{BayesMfp} object
  with the same attributes, which are

  \item{numVisited}{the number of models that have been visited
  (exhaustive search) or cached (model sampling)}
  \item{inclusionProbs}{BMA inclusion probabilities for all uncertain
    covariates}
  \item{linearInclusionProbs}{BMA probabilities for exactly linear
  inclusion of FP covariates} 
  \item{logNormConst}{the (estimated) log normalizing constant \eqn{f
      (D)}}
  \item{chainlength}{length of the Markov chain, only present if \code{method = "sampling"}}
  \item{call}{the original call}
  \item{formula}{the formula by which the appropriate untransformed
    design matrix can be extracted}
  \item{x}{the shifted and scaled design matrix for the data}
  \item{xCentered}{the column-wise centered x}
  \item{y}{the response vector}
  \item{yMean}{the mean of the response values}
  \item{SST}{sum of squares total}
  \item{indices}{a list with components that describe the positions of
    uncertain covariate groups, fractional polynomial terms and fixed
    variables in the design matrix}
  \item{termNames}{a list of character vectors containing the names of
    uncertain covariate groups, fractional polynomial terms and fixed
    variables}
  \item{shiftScaleMax}{matrix with 4 columns containing preliminary
    transformation parameters, maximum degrees and cardinalities of the
    powersets of the fractional polynomial terms}
  \item{priorSpecs}{the utilized prior specifications}
  \item{randomSeed}{if a seed existed at function call
  (\code{get(".Random.seed", .GlobalEnv)}), it is saved here} 
}

\note{\code{logNormConst} may be unusable due to necessary conversion
  from long double to double!

  Various methods for posterior summaries are available.
}

\seealso{
  \link{BayesMfp Methods}, \code{\link{BmaSamples}}
}

\examples{
## generate some data
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

## run an exhaustive model space evaluation with a flat model prior and
## a uniform prior (a = 4) on the shrinkage factor t = g/(1 + g):
test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")
test

## now the same with a *dependent* model prior:
test2 <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
		 priorSpecs = list(a = 4, modelPrior = "dependent"),
                 method="exhaustive")
test2
}

\keyword{regression}
