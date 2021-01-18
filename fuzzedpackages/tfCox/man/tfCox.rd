\name{tfCox}
\alias{tfCox}
\title{Fit the additive trend filtering Cox model with a range of tuning parameters}
\description{Fit additive trend filtering Cox model where each component function is estimated to be piecewise constant or polynomial.
}
\usage{
tfCox(dat, ord=0, alpha=1, lambda.seq=NULL, discrete=NULL, n.lambda=30,
lambda.min.ratio = 0.01, tol=1e-6, niter=1000, stepSize=25, backtracking=0)
}
\arguments{
  \item{dat}{
A list that contains \code{time}, \code{status} and \code{X}. \code{time} is failure or censoring time, \code{status} is failure indicator with 1 indicating failure and 0 indicating censoring, and \code{X} is n x p design matrix and may have p > n. Missing data are not allowed in \code{time}, \code{status} and \code{X}. \code{X} should be numeric.
}
  \item{ord}{
The polynomial order of the trend filtering fit; a non-negative interger (\code{ord>= 3} is not recommended). For instance, \code{ord=0} will produce piewise constant fit, \code{ord=1} will produce piewise linear fit, and \code{ord=2} will produce piewise quadratic fit.
}
  \item{alpha}{
The trade-off between trend filtering penalty and group lasso penalty. It must be in [0,1]. \code{alpha=1} corresponds to the case with only trend filtering penalty to produce piecewise polynomial, and \code{alpha=0} corresponds to the case with only group lasso penalty to produce sparsity of the functions. \code{alpha} between 0 and 1 is the tradeoff between the strength of these two penalties. For p < n, we suggest using 1.
}
  \item{lambda.seq}{
A vector of non-negative tuning parameters. If provided, \code{lambda.seq} should be a decreasing sequence of values since \code{tfCox} uses warm starts for speed. If \code{lambda.seq=NULL}, the default will calculate \code{lambda.seq} using \code{lambda.min.ratio} and \code{n.lambda}.
}
  \item{discrete}{
A vector of covariate/feature indice that are discrete. Discrete covariates are not penalized in the model. Default \code{NULL} means that none of the covariates are discrete thus all covariates will be penalized in the model.
}
  \item{n.lambda}{
The number of lambda values to consider and the default is 30.
}
  \item{lambda.min.ratio}{
Smallest value for lambda.seq, as a fraction of the maximum lambda value, which is the smallest value such that the penalty term is zero. The default is 0.01.
}
  \item{tol}{
Convergence criterion for estimates.
}
  \item{niter}{
Maximum number of iterations.
}
  \item{stepSize}{
Initial step size. Default is 25.
}
    \item{backtracking}{
Whether backtracking should be used 1 (TRUE) or 0 (FALSE). Default is 0 (FALSE). 
}
}

\value{
An object with S3 class "tfCox".
\item{ord}{the polynomial order of the trend filtering fit. Specified by user (or default).  }
\item{alpha}{as specified by user (or default). }
\item{lambda.seq}{vector of lambda values considered. }
\item{theta.list}{list of estimated theta matrices of dimension n x p. Each component in the list corresponds to the fit from \code{lambda.seq}. }
\item{num.knots}{vector of number of knots of the estimated theta. Each component corresponds to the fit from \code{lambda.seq}. }
\item{num.nonsparse}{vector of proportion of non-sparse/non-zero covariates/features. Each component corresponds to the fit from \code{lambda.seq}. }
\item{dat}{as specified by user. }
}

\details{
The optimization problem has the form

\deqn{ l(\theta)+\alpha\lambda\sum_{j=1}^p |D_jP_j\theta_j|_1+(1-\alpha)\lambda\sum_{j=1}^p|\theta_j|_2 }

where \eqn{l} is the loss function defined as the negative log partial likelihood divided by n, and \eqn{\alpha} provides a trade-off between trend filtering penalty and group lasso penalty. Covariate matrix \code{X} is not standardized before solving the optimization problem.

}

\author{
Jiacheng Wu
}

\references{
Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758
}

\seealso{
\code{\link{summary.tfCox}}, \code{\link{predict.tfCox}}, \code{\link{plot.tfCox}}, \code{\link{cv_tfCox}}
}

\examples{
###################################################################
#constant trend filtering (fused lasso) with adaptively chosen knots
#generate data from simulation scenario 1 with piecewise constant functions
set.seed(1234)
dat = sim_dat(n=100, zerof=0, scenario=1)

#fit piecewise constant for alpha=1 and a range of lambda
fit = tfCox(dat, ord=0, alpha=1)
summary(fit)
#plot the fit of lambda index 15 and the first predictor
plot(fit, which.lambda=15, which.predictor=1)

#cross-validation to choose the tuning parameter lambda with fixed alpha=1
cv = cv_tfCox(dat, ord=0, alpha=1, n.fold=2)
summary(cv)
cv$best.lambda
#plot the cross-validation curve
plot(cv)

#fit the model with the best tuning parameter chosen by cross-validation
one.fit = tfCox(dat, ord=0, alpha=1, lambda.seq=cv$best.lambda)
#predict theta from the fitted tfCox object
theta_hat = predict(one.fit, newX=dat$X, which.lambda=1)

#plot the fitted theta_hat (line) with the true theta (dot)
for (i in 1:4) {
  ordi = order(dat$X[,i])
  plot(dat$X[ordi,i], dat$true_theta[ordi,i],
    xlab=paste("predictor",i), ylab="theta" )
  lines(dat$X[ordi,i], theta_hat[ordi,i], type="s")
}

\donttest{
#################################################################
#linear trend filtering with adaptively chosen knots
#generate data from simulation scenario 3 with piecewise linear functions
set.seed(1234)
dat = sim_dat(n=100, zerof=0, scenario=3)

#fit piecewise constant for alpha=1 and a range of lambda
fit = tfCox(dat, ord=1, alpha=1)
summary(fit)
#plot the fit of lambda index 15 and the first predictor
plot(fit, which.lambda=15, which.predictor=1)

#cross-validation to choose the tuning parameter lambda with fixed alpha=1
cv = cv_tfCox(dat, ord=1, alpha=1, n.fold=2)
summary(cv)
#plot the cross-validation curve
plot(cv)

#fit the model with the best tuning parameter chosen by cross-validation
one.fit = tfCox(dat, ord=1, alpha=1, lambda.seq=cv$best.lambda)
#predict theta from the fitted tfCox object
theta_hat = predict(one.fit, newX=dat$X, which.lambda=1)

#plot the fitted theta_hat (line) with the true theta (dot)
for (i in 1:4) {
  ordi = order(dat$X[,i])
  plot(dat$X[ordi,i], dat$true_theta[ordi,i],
       xlab=paste("predictor",i), ylab="theta" )
  lines(dat$X[ordi,i], theta_hat[ordi,i], type="l")
}
}

}
