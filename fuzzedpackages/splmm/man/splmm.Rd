\name{splmm}
\alias{splmm}
\alias{splmm.default}
\title{Function to fit linear mixed-effects model with double penalty for fixed effects and random effects}
\usage{
splmm(x, y, z, grp, lam1, lam2, nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
      penalty.L=c("lasso","scad"),CovOpt=c("nlminb","optimize"),
      standardize=TRUE,control=splmmControl())

\method{splmm}{default}(x, y, z, grp, lam1, lam2, nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
      penalty.L=c("lasso","scad"),CovOpt=c("nlminb","optimize"),
      standardize=TRUE,control=splmmControl())}

\arguments{
  \item{x}{matrix of dimension N x p including the fixed-effects covariables. An intercept has to be included in the first column as (1,...,1).}
  \item{y}{response variable of length N.}
  \item{z}{random effects matrix of dimension N x q. It has to be a matrix, even if q=1.}
  \item{grp}{grouping variable of length N}
  \item{lam1}{regularization parameter for fixed effects penalization.}
  \item{lam2}{regularization parameter for random effects penalization.}
  \item{nonpen.b}{Index of indices of fixed effects not penalized. The default value is 1, which means the fixed intercept is not penalized}
  \item{nonpen.L}{Index of indices of random effects not penalized. The default value is 1, which means the random intercept is not penalized}
  \item{penalty.b}{The penalty method for fixed effects penalization. Currently available options include LASSO penalty and SCAD penalty. }
  \item{penalty.L}{The penalty method for fixed effects penalization. Currently available options include LASSO penalty and SCAD penalty. }
  \item{CovOpt}{which optimization routine should be used for updating the variance parameter. The available options include optimize and nlminb. nlminb uses the estimate of the last iteration as a starting value. nlminb is faster if there are many Gauss-Seidel iterations.}
  \item{standardize}{A logical parameter specifying whether the fixed effects matrix x and random effects matrix z should  be standardized such that each column has mean 0 and standard deviation 1. The default value is \code{TRUE}}
  \item{control}{control parameters for the algorithm and the Armijo Rule, see \code{'splmmControl'} for the details}
}

\description{
All the details of the algorithm can be found in the manuscript.
}


\value{A \code{'splmm'} object is returned, for which
  \code{coef},\code{resid}, \code{fitted},
  \code{print}, \code{summary} methods exist.


  \item{data}{data set used for fitting the model, as a list with four components: x, y, z, grp (see
    above)}
  \item{coefInit}{list of the starting values for beta, random effects covariance structure, and variance structure}
  \item{penalty.b}{The penalty method for fixed effects penalization.}
  \item{penalty.L}{The penalty method for random effects penalization.}
  \item{nonpen.b}{Index of indices of fixed effects not penalized.}
  \item{nonpen.L}{Index of indices of random effects not penalized.}
  \item{lambda1}{regularization parameter for fixed effects penalization scaled by the number of subjects.}
  \item{lambda2}{regularization parameter for random effects penalization the number of subjects.}
  \item{sigma}{standard deviation \eqn{\hat{\sigma}} of the errors}
  \item{D}{The estimates of the random effects covariance matrix \eqn{\hat{D}}. }
  \item{Lvec}{Vectorized \eqn{\hat{L}}, the lower triangular matrix of \eqn{\hat{D}} from Cholesky Decomposition. }
  \item{coefficients}{estimated fixed-effects coefficients \eqn{\hat{\beta}}}
  \item{random}{vector with random effects, sorted by groups}
  \item{ranef}{vector with random effects, sorted by effect}
  \item{u}{vector with the standardized random effects, sorted by effect}
  \item{fixef}{estimated fixed-effects coeffidients \eqn{\hat{\beta}}}
  \item{fitted.values}{The fitted values \eqn{\hat{y} = \hat{X} \beta +
  Z \hat{b}_i}}
  \item{residuals}{raw residuals \eqn{y-\hat{y}}}
  \item{corD}{Correlation matrix of the random effects}
  \item{logLik}{value of the log-likelihood function}
  \item{deviance}{deviance=-2*logLik}
  \item{npar}{Number of parameters. Corresponds to the cardinality
  of the set of nonzero \code{coefficients} plus the number of nonzero variance in \code{D}}
  \item{aic}{AIC}
  \item{bic}{BIC}
  \item{bbic}{Modified BIC defined by Wang et al (2009)}
  \item{ebic}{Extended BIC defined by Chen and Chen (2008)}
  \item{converged}{Does the algorithm converge? 0: correct convergence ;
  an odd number means that maxIter was reached ; an even number means
  that the Armijo step was not succesful. For each unsuccessfull Armijo
  step, 2 is added to converged. If converged is large compared to the
  number of iterations \code{counter}, you may increase maxArmijo.}
  \item{counter}{The number of iterations used.}
  \item{stopped}{logical indicating whether the algorithm stopped due to too many parameters, if yes need to increase \code{lam1} or \code{lam2}}
  \item{CovOpt}{optimization routine}
  \item{control}{see \code{splmmControl}}
  \item{objective}{Value of the objective function at the final estimates}

  \item{call}{call}


}
\examples{
\donttest{

## Use splmm on the Kenya school cognitive data set


data(cognitive)

x <- model.matrix(ravens ~schoolid+treatment+year+sex+age_at_time0
                  +height+weight+head_circ+ses+mom_read+mom_write
                  +mom_edu, cognitive)
z <- x

fit <- splmm(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1=0.1,
lam2=0.1,penalty.b="lasso", penalty.L="lasso")
summary(fit)
}

}


