\name{gen.data}
\alias{gen.data}
\title{Generate simulated data}
\description{
   Generate data for simulations under the generalized linear model and Cox model.
}
\usage{
  gen.data(n, p, family, K, rho = 0, sigma = 1, beta = NULL, censoring = TRUE,
           c = 1, scal)
}
\arguments{
  \item{n}{The number of observations.}
  \item{p}{The number of predictors of interest.}
   \item{family}{The distribution of the simulated data. "\code{gaussian}" for gaussian data."\code{binomial}" for binary data. "\code{cox}" for survival data }
  \item{K}{The number of nonzero coefficients in the underlying regression model.}
  \item{rho}{A parameter used to characterize the pairwise correlation in predictors. Default is 0.}
  \item{sigma}{A parameter used to control the signal-to-noise ratio. For linear regression, it is the error variance \eqn{\sigma^2}. For logistic regression and Cox's model, the larger the value of sigma, the higher the signal-to-noise ratio. }
  \item{beta}{The coefficient values in the underlying regression model.}
  \item{censoring}{Whether data is censored or not. Default is TRUE}
  \item{c}{The censoring rate. Default is 1.}
  \item{scal}{A parameter in generating survival time based on the Weibull distribution. Only used for the "\code{cox}" family. }
}
\value{
A list with the following components: x, y, Tbeta.
 \item{x}{Design matrix of predictors.}
 \item{y}{Response variable}
  \item{Tbeta}{The coefficients used in the underlying regression model.}
}
\details{
 For the design matrix \eqn{X}, we first generate an n x p random Gaussian matrix \eqn{\bar{X}} whose entries are i.i.d. \eqn{\sim N(0,1)} and then normalize its columns to the \eqn{\sqrt n} length. Then the design matrix \eqn{X} is generated with \eqn{X_j = \bar{X}_j + \rho(\bar{X}_{j+1}+\bar{X}_{j-1})} for \eqn{j=2,\dots,p-1}.

 For "\code{gaussian}" family, the data model is
  \deqn{Y = X \beta + \epsilon,  where  \epsilon \sim N(0, \sigma^2 ).}
 The underlying regression coefficient \eqn{\beta} has uniform distribution [m, 100m], \eqn{m=5 \sqrt{2log(p)/n}.}

 For "\code{binomial}" family, the data model is
  \deqn{Prob(Y = 1) = exp(X \beta)/(1 + exp(X \beta))}
 The underlying regression coefficient \eqn{\beta} has uniform distribution [2m, 10m], \eqn{m = 5\sigma \sqrt{2log(p)/n}.}

 For "\code{cox}" family, the data model is
 \deqn{T = (-log(S(t))/exp(X \beta))^(1/scal),}
 The centerning time \code{C} is generated from uniform distribution [0, c], then we define the censor status as \eqn{\delta = I{T <= C}, R = min{T, C}}.
  The underlying regression coefficient \eqn{\beta} has uniform distribution [2m, 10m], \eqn{m = 5\sigma \sqrt{2log(p)/n}.}
}

\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S. and Wang, X. (2020). BeSS: An R Package for Best Subset Selection in Linear, Logistic and Cox Proportional Hazards Models, \emph{Journal of Statistical Software}, Vol. 94(4).  doi:10.18637/jss.v094.i04.
}
\examples{

# Generate simulated data
n <- 500
p <- 20
K <-10
sigma <- 1
rho <- 0.2
data <- gen.data(n, p, family = "gaussian", K, rho, sigma)

# Best subset selection
fit <- bess(data$x, data$y, family = "gaussian")


}
