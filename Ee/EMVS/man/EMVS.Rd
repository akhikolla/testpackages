\name{EMVS}
\alias{EMVS}
\title{Bayesian Variable Selection using EM Algorithm}
\description{
EMVS is a fast deterministic approach to identifying sparse high posterior models for Bayesian variable selection under spike-and-slab priors in linear regression.  
EMVS performs dynamic posterior exploration, which outputs a solution path computed at a grid
  of values for the spike variance parameter \code{v0}. 
  }
\usage{
EMVS(Y, X, v0, v1, type = c("betabinomial", "fixed"), independent = TRUE,  
      beta_init, sigma_init, epsilon = 10^(-5), temperature, theta, a, b, v1_g,
    direction=c("backward", "forward", "null"), standardize = TRUE, log_v0 = FALSE)
}
\arguments{
  \item{Y}{Vector of continuous responses (n x 1). The  responses are expected to be centered.}
  \item{X}{Matrix of regressors (n x p). Continous predictors are expected to be standardized to have mean zero and standard deviation one.
}
  \item{v0}{Spike variance parameters. Either a numeric value for a single run or a sequence of increasing values for dynamic posterior exploration.}
  \item{v1}{Slab variance parameter. Needs to be greater than v0.}
  \item{type}{Type of the prior distribution over the model space: \code{type="betabinomial"} for the betabinomial prior with shape parameters \code{a} and \code{b}, \code{type="fixed"} for the Bernoulli prior with a fixed inclusion probability \code{theta}.}
  \item{independent}{If TRUE, the regression coefficients and the error variance are taken to be independent a priori (default). If FALSE, a conjugate prior is used as in Rockova and George (2014).}
  \item{beta_init}{Vector (p x 1) of initial values for the regression parameters \code{beta}. If missing, a default vector of starting values obtained as a limiting case of deterministic annealing used
\deqn{beta^0=[X'X+0.5(1/v1+1/v0)I_p]^{-1}X'Y.}}
  \item{sigma_init}{Initial value for the residual variance parameter.}
  \item{epsilon}{Convergence margin parameter.  The computation at each \code{v0} is terminated when 
      \deqn{ ||beta^{k+1}-beta^{k}||_2<epsilon.}}
  \item{temperature}{Temperature parameter for deterministic annealing. If missing, a default value \code{temperature=1} used.}
  \item{theta}{Prior inclusion probability for \code{type="fixed"}.}
  \item{a,b}{Scale parameters of the beta distribution for \code{type="betabinomial"}.}
  \item{v1_g}{Slab variance parameter value for the g-function. If missing, a default value \code{v1} used.}
  \item{direction}{Direction of the sequential reinitialization in dynamic posterior exploration. The default is \code{direction="backward"} - this initializes the first computation at \code{beta_init} using the largest value of \code{v0} and uses the resulting output as a warm start for the next largest value \code{v0} in a backward direction (i.e. from the largest to the smallest value of \code{v0}). The option \code{direction="forward"} proceeds from the smallest value of \code{v0} to the largest value of \code{v0}, using the output from the previous solution as a warm start for the next. \code{direction = "null"} re-initializes at \code{beta_init} for each \code{v0}.}
  \item{standardize}{If \code{TRUE} (default), the design matrix X is standardized (mean zero and variance \code{n}). }
  \item{log_v0}{If TRUE, the v0s are displayed on the log scale in \code{EMVSplot}.}
 }

\details{
	An EM algorithm is applied to find posterior modes of the regression parameters in linear models under spike and slab priors.
Variable selection is performed by threshodling the posterior modes to obtain models \code{gamma} with high posterior probability \code{P(gamma|Y)} . The spike variance \code{v0} can be altered to obtain models with various degrees of sparsity. The slab variance is set to a fixed value \code{v1>v0}. The thresholding is based on the conditional posterior probabilities of inclusion, which are outputed of the procedure. Variables are included as long as their inclusion probability is above  \code{0.5}.
Dynamic exploration is achieved by considering a sequence of increasing spike variance parameters \code{v0}. For each v0, a candidate model is obtained.
For the conjugate prior case, the best model is then picked according to a criterion ("log g-function"), which equals to the log of the posterior model probability  up to a constant
\deqn{
  log g(gamma)=log P(gamma| Y) + C.
} 
Independent and sequential initializations are implemented. Sequential initialization uses previously found modes are warm starts in both forward and backward direction of the given sequence of v0 values.}

\value{
A list object, for which \code{EMVSplot} and \code{EMVSbest} functions exist.

  \item{betas}{Matrix of estimated regression coefficients (posterior modal estimates) of dimension  (L x p), where L is the length of v0.}
  \item{log_g_function}{Vector (L x 1) of log posterior model probabilities (up to a constant) of subsets found for each v0. (Only available for \code{independent = FALSE}).}
  \item{intersects}{Vector (L x 1) of posterior weighted intersection points between spike and slab components.}
  \item{sigmas}{Vector (L x 1) of estimated residual variances.}
  \item{v1}{Slab variance parameter values used.}
  \item{v0}{Spike variance parameter values used.}
  \item{niters}{Vector (L x 1) of numbers of interations until convergence for each v0}
  \item{prob_inclusion}{A matrix (L x p) of conditional inclusion probabilities. Each row corresponds to a single v0 value. }
  \item{type}{Type of the model prior used.}
  \item{type}{Type of initialization used, \code{type="null"} stands for the default cold start.}
  \item{theta}{Vector (L x 1) of estimated inclusion probabilities for \code{type="betabinomial"}.}
}
\author{Veronika Rockova\cr 
Maintainer: Veronika Rockova \email{veronika.rockova@chicagobooth.edu}}
\references{Rockova, V. and George, E. I.  (2014)
  \emph{EMVS: The EM Approach to Bayesian Variable Selection},   
\url{http://amstat.tandfonline.com/doi/abs/10.1080/01621459.2013.869223?journalCode=uasa20#preview}\cr
  \emph{Journal of the American Statistical Association}
}


\seealso{\code{EMVSplot}, \code{EMVSsummary}, \code{EMVSbest}}

\examples{
# Linear regression with p>n variables
library(EMVS)

n = 100
p = 1000
X = matrix(rnorm(n * p), n, p)
beta = c(1.5, 2, 2.5, rep(0, p-3))
Y = X[,1] * beta[1] + X[,2] * beta[2] + X[,3] * beta[3] + rnorm(n)

# conjugate prior on regression coefficients and variance
v0 = seq(0.1, 2, length.out = 20)
v1 = 1000
beta_init = rep(1, p)
sigma_init = 1
a = b = 1
epsilon = 10^{-5}

result = EMVS(Y, X, v0 = v0, v1 = v1, type = "betabinomial", 
independent = FALSE, beta_init = beta_init, sigma_init = sigma_init,
epsilon = epsilon, a = a, b = b)

EMVSplot(result, "both", FALSE)

EMVSbest(result)

# independent prior on regression coefficients and variance
v0 = exp(seq(-10, -1, length.out = 20))
v1 = 1
beta_init = rep(1,p)
sigma_init = 1
a = b = 1
epsilon = 10^{-5}

result = EMVS(Y, X, v0 = v0, v1 = v1, type = "betabinomial", 
independent = TRUE, beta_init = beta_init, sigma_init = sigma_init,
epsilon = epsilon, a = a, b = b, log_v0 = TRUE)

EMVSplot(result, "both", FALSE)

EMVSbest(result)
}
\keyword{Bayesian variable selection}
\keyword{Spike and slab}
