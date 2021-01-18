\name{fdaLm}
\alias{fdaLm}
\title{Linear mixed-effects model for functional data}
\description{Fits variance and smoothing parameters, and possibly also
  Box-Cox transformation, by maximum restricted likelihood. Estimate
  fixed parameters, predict random effects, and predict serial
  correlated effect at point of maximum restricted likelihood. Linear
  models for fixed and random effects may be global or marginal over
  sample times.}
\usage{
fdaLm(formula, data, design, boxcox = NULL, G = 1, lambda = 1, nlSearch
= TRUE, K.order = 1, D.order = NULL, Fleft = "tied", Fright = "tied",
left = NULL, right = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    A multiple formula of the type \code{Y|id ~ fixed|random}. Here
    \code{Y} is the response variable, \code{id} is a factor 
    separating the samples, \code{fixed} is a linear model for the fixed
    effect, and \code{random} is a linear model for the random
    effect.}
  \item{data}{
    An optional data frame containing the variables. See details below.}
  \item{design}{
    An optional data frame containing the design variables in the
    specification of the fixed and the random effects. See details below.}
  \item{boxcox}{
    The power parameter in the scale invariant Box-Cox transformation.
    If \code{NULL} (default), then no transformation is performed. If a
    numeric value is provided, then a scale invariant Box-Cox
    transformation of the response variable is performed. The numeric
    value is either used as it is (\code{nlSearch=FALSE}) or as the
    starting point for a non-linear optimization (\code{nlSearch=TRUE}.)}
  \item{G}{
    Variance of the random effects. Present implementation only allows
    for independent random effects, i.e. \code{G} is scalar. Used
    depending on \code{nlSearch} as described above.}
  \item{lambda}{
    Start value for the \code{lambda} parameter describing the
    L-operator. Presently the following forms are implemented: If
    \code{K.order} is odd, then \code{lambda} may have length=1
    corresponding to \code{L=-lambda[1]*D^(2*K.order)}, or length=2
    corresponding to \code{L=-lambda[1]*D^(2*K.order)+lambda[2]}. If
    \code{K.order} is even, then \code{lambda} may have length=1
    corresponding to \code{L=-lambda[1]*D^(2*K.order)}, length=2
    corresponding to
    \code{L=-lambda[1]*D^(2*K.order)+lambda[2]*D^K.order}, or length=3
    corresponding to
    \code{L=-lambda[1]*D^(2*K.order)+lambda[2]*D^K.order+lambda[3]}.
    Used depending on \code{nlSearch} as described above.
    All coefficients must be non-negative, and the leading coefficient
    \code{lambda[1]} must be strictly positive. Coefficients equal to
    zero are kept fixed at zero in the non-linear optimization.}
  \item{nlSearch}{
    If \code{TRUE} (default), then a non-linear optimization of the
    parameters \code{boxcox}, \code{G}, \code{lambda} is performed
    (present implementation uses \code{nlminb}). If \code{FALSE},
    then the initial values of the non-linear parameters are used.}
  \item{K.order}{
    The order of the K-operator.}
  \item{D.order}{
    The requested order of derivatives of the prediction of the serial
    correlated effect \code{xBLUP}. If \code{NULL} (default), then
    \code{D.order} is set to the maximal recommended order
    \code{K.order}.}
  \item{Fleft}{
    Specification of the \code{K.order} boundary conditions at the left
    limit of the sampling interval. Value \code{"tied"} (default) gives
    bridge-type conditions. Value \code{"open"} shifts up the
    bridge-type conditions one differential order, hence removing the
    restriction on the level (corresponding to the open end of a
    Brownian motion). Otherwise arbitrary linear boundary conditions may
    be specified as a matrix with dimension
    (\code{K.order},2*\code{K.order}).}
  \item{Fright}{
    Similarly for the \code{K.order} boundary conditions at the right
    limit of the sampling interval.}
  \item{left}{
    Left limit of the sampling interval. If \code{NULL} (default), then
    \code{left} is set to 0.}
  \item{right}{
    Right limit of the sampling interval. If \code{NULL} (default), then
    \code{right} is set to the number of sampling points. Thus, the
    default values of \code{left} and \code{right} give sampling distance
    equal to 1.}
%  \item{lag.max}{
%    The maximal number of lags used in the computation of the
%    theoretical counterpart of the empirical auto covariance function.}
}
\details{
  The response variable \code{Y} is taken from the data frame
  \code{data} (subsidiary the parent environment). If there is more
  than one sample, then the responses must be stacked sample-wise on top
  of each other. The sample identifier \code{id} is sought for in both
  data frames \code{data} and \code{design} (subsidiary the parent
  environment). The primarily function of the identifier is to decide
  the number of samples. But if \code{id} is present in both data
  frames, and if there is more that one sample, then this variable is
  also used to match the reponse vector to the design variables
  (i.e. these need not appear in the same order).

  The design variables \code{fixed} and \code{random} for the fixed and
  the random effects are taken from the data frame \code{design}
  (subsidiary the parent environment), subsidiary from the data frame
  \code{data} (subsidiary the parent environment).

  If the number of observations in the design variables equal the total
  number of response observations, then a global ANOVA is performed. If the
  number of observations in the design variables equal the number of
  sample points, then a marginal ANOVA is performed.
}
\note{
  If the real value of the left most eigenvalues are non-positive, and
  if the real value of the right most eigenvalues are non-negative,
  then the underlying algorithm is numerical stable. This will always be
  the situation for the present restriction of the L-operator. 

  If \code{lambda} has length=1, then it may also be interpreted as the
  smoothing parameter in the penalized likelihood framework.

  If \code{D.order} is chosen larger than \code{K.order}, this number of
  derivaties are also computed during the non-linear optimization. This
  might slow down the computation speed a little bit.
  }
\value{
  A list with components
  \item{logLik}{Minus twice the log restricted likelihood taken at the
    estimates.}
  \item{ANOVA}{Specifies whether fixed and random effects were estimated
    globally (\code{global}) or marginally (\code{marginal}).}
  \item{nlSearch}{Specifies whether non-linear optimization was
    performed (\code{TRUE} / \code{FALSE}).}
  \item{counts}{Number of computations of the negative log likelihood.}
  \item{boxcoxHat}{Maximum restricted likelihood estimate for the power
    parameter in the scale invariant Box-Cox transformation. Equal to
    \code{not done} if the Box-Cox transformation is not used.}
  \item{Ghat}{Maximum restricted likelihood estimate for the variance
    matrix of the random effects.}
  \item{lambdaHat}{Maximum restricted likelihood estimate for the lambda
    parameter describing the L-operator.}
  \item{sigma2hat}{Maximum restricted likelihood estimate for the noise
    variance.}
  \item{betaHat}{For global ANOVA a vector with estimate for the fixed
    effect. For marginal ANOVA a matrix with estimate for the fixed
    effects.}
  \item{uBLUP}{For global ANOVA a vector with prediction of the random
    effect. For marginal ANOVA a matrix with predictions of the random
    effects.}
  \item{xBLUP}{Array with predictions of serial correlated effects. The
    dimension is (sample length,sample numbers,1+\code{D.order}).}
  \item{condRes}{Matrix of conditional residuals. The dimension is
    (sample length,sample numbers).}
  \item{betaVar}{Variance matrix of fixed effect estimate.}
%  \item{uVar}{Variance matrix of random effect prediction.}
}
\author{Bo Markussen <bomar@math.ku.dk>}
\seealso{See also \code{\link{findRoots}} and \code{\link{dataTrans}}.}
\examples{
# ---------------------
# Using a fixed effect
# ---------------------
x <- seq(0,2*pi,length.out=200)
y.true <- sin(x)+x
y.obs <- y.true + rnorm(200)
est0 <- fdaLm(y.obs~0,Fright="open",right=2*pi)
est1 <- fdaLm(y.obs~0+x,Fright="open",right=2*pi)
plot(x,y.obs,main="Estimating the sum of a line and a curve")
lines(x,y.true,lty=2)
lines(x,est0$xBLUP[,1,1],col=2)
lines(x,est1$betaHat*x+est1$xBLUP[,1,1],col=3)
legend("topleft",c("True curve","Smooth","Line + smooth"),col=1:3,lty=c(2,1,1))

# --------------------------
# Including a random effect
# --------------------------
# Build data frame
test.frame <- data.frame(y=rnorm(50),sample=factor(rep(1:5,each=10)),
                         x=rep(0:9,times=5),
                         f=factor(rnorm(50) < 0,labels=c("a","b")),
                         j=factor(rnorm(50) < 0,labels=c("A","B")))
test.frame$y <- test.frame$y + 2 +
    3*(test.frame$f=="a")*test.frame$x + 5*(test.frame$f=="b")*test.frame$x +
(-10)*(test.frame$j=="A") + 10*(test.frame$j=="B")
# This is the model 'y|sample ~ f:x|j' with intercept=2, slopes (3,5),
# and random effects (-10,10)
est <- fdaLm(y|sample ~ f:x|0+j,data=test.frame)
print(est)
}
\keyword{ model }
\keyword{ inference }
