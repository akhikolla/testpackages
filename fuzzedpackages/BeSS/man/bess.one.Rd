\name{bess.one}
\alias{bess.one}
\title{Best subset selection with a specified model size}
\description{
  Best subset selection with a specified model size for generalized linear models and Cox's proportional hazard model.
}
\usage{
bess.one(x, y, family = c("gaussian", "binomial", "cox"),
         s = 1,
         max.steps = 15,
         glm.max = 1e6,
         cox.max = 20,
         factor = NULL,
         weights = rep(1,nrow(x)),
         normalize = TRUE)
}
\arguments{
  \item{x}{Input matrix,of dimension n x p; each row is an observation vector.}
  \item{y}{Response variable, of length n. For family = "\code{gaussian}", \code{y} should be a vector with continuous values. For family = "\code{binomial}", \code{y} should be a factor with two levels. For family = "\code{cox}", \code{y} should be a two-column matrix with columns named 'time' and 'status'.}
  \item{s}{Size of the selected model.It controls number of nonzero coefiicients to be allowed in the model.}
  \item{family}{One of the ditribution function for GLM or Cox models. Either "\code{gaussian}", "\code{binomial}", or "\code{cox}", depending on the response.}
  \item{max.steps}{The maximum number of iterations in the primal dual active set algorithm. In most cases, only a few steps can gurantee the convergence. Default is 15.}
  \item{glm.max}{The maximum number of iterations for solving the maximum likelihood problem on the active set. It occurs at each step in the primal dual active set algorithm. Only used in the logistic regression for family = "\code{binomial}". Default is \eqn{1e+6}.}
  \item{cox.max}{The maximum number of iterations for solving the maximum partial likelihood problem on the active set. It occurs at each step in the primal dual active set algorithm. Only used in Cox model for family  = "\code{cox}". Default is 20.}
  \item{weights}{Observation weights. Default is 1 for each observation}
  \item{factor}{Which variable to be factored. Should be NULL or a numeric vector.}
  \item{normalize}{Whether to normalize \code{x} or not. Default is TRUE.}
}
\value{
A list with class attribute '\code{bess.one}' and named components:
  \item{type}{Types of the model: "\code{bess_gaussian}" for linear model, "\code{bess_binomial}" for logistic model and "\code{bess_cox}" for Cox model}
  \item{beta}{The best fitting coefficients with the smallest loss function given the model size \code{s}.}
  \item{lambda}{The estimated lambda value in the Lagrangian form of the best subset selection problem with model size \code{s}.}
  \item{bestmodel}{The best fitted model, the class of which is "lm", "glm" or "coxph"}
  \item{deviance}{The value of \eqn{-2*logL(\beta)}.}
  \item{nulldeviance}{The value of \eqn{-2*logL(\beta)} for null model.}
  \item{factor}{Which variable to be factored. Should be NULL or a numeric vector.}
}

\details{
  Given a model size \eqn{s}, we consider the following best subset selection problem:
  \deqn{\min_\beta -2 logL(\beta) ;{ s.t.} \|\beta\|_0 = s.}
  In the GLM case, \eqn{logL(\beta)} is the log-likelihood function; In the Cox model, \eqn{logL(\beta)} is the log parital likelihood function.

  The best subset selection problem is solved by the primal dual active set algorithm, see Wen et al. (2017) for details. This algorithm utilizes an active set updating strategy via primal and dual vairables and fits the sub-model by exploiting the fact that their support set are non-overlap and complementary.
}

\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S. and Wang, X. (2020). BeSS: An R Package for Best Subset Selection in Linear, Logistic and Cox Proportional Hazards Models, \emph{Journal of Statistical Software}, Vol. 94(4).  doi:10.18637/jss.v094.i04.
}


\seealso{
\code{\link{bess}}, \code{\link{plot.bess}},
\code{\link{predict.bess}}.
}


\examples{

#--------------linear model--------------#
# Generate simulated data

n <- 500
p <- 20
K <-10
sigma <- 1
rho <- 0.2
data <- gen.data(n, p, family = "gaussian", K, rho, sigma)


# Best subset selection
fit1 <- bess.one(data$x, data$y, s = 10, family = "gaussian", normalize = TRUE)
#coef(fit1,sparse=TRUE)
bestmodel <- fit1$bestmodel
#summary(bestmodel)

## Not run: 
#--------------logistic model--------------#

# Generate simulated data
data <- gen.data(n, p, family = "binomial", K, rho, sigma)

# Best subset selection
fit2 <- bess.one(data$x, data$y, family = "binomial", s = 10, normalize = TRUE)
bestmodel <- fit2$bestmodel
#summary(bestmodel)

#--------------cox model--------------#

# Generate simulated data
data <- gen.data(n, p, K, rho, sigma, c=10, family="cox", scal=10)

# Best subset selection
fit3 <- bess.one(data$x, data$y, s = 10, family = "cox", normalize = TRUE)
bestmodel <- fit3$bestmodel
#summary(bestmodel)

#----------------------High dimensional linear models--------------------#

p <- 1000
data <- gen.data(n, p, family = "gaussian", K, rho, sigma)

# Best subset selection
fit <- bess.one(data$x, data$y, s=10, family = "gaussian", normalize = TRUE)

#---------------prostate---------------#
data("prostate")
x = prostate[,-9]
y = prostate[,9]

fit.ungroup = bess.one(x, y, s=5)
fit.group = bess.one(x, y, s=5, factor = c("gleason"))

#---------------SAheart---------------#
data(SAheart)
y = SAheart[,5]
x = SAheart[,-5]
x$ldl[x$ldl<5] = 1
x$ldl[x$ldl>=5&x$ldl<10] = 2
x$ldl[x$ldl>=10] = 3

fit.ungroup = bess.one(x, y, s=5, family = "binomial")
fit.group = bess.one(x, y, s=5, factor = c("ldl"), family = "binomial")
## End(Not run)
}
