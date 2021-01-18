\name{tfCox_choose_lambda}
\alias{tfCox_choose_lambda}
\title{
Choose the tuning parameter lambda using training and testing dataset
}
\description{
Fit additive trend filtering Cox model where each component function is estimated to be piecewise constant or polynomial. Tuning parameter is selected via training and testing dataset described in Wu and Witten (2019). Training data is used to build the model, and testing data is used for selecting tuning parameter based on log likelihood. It is a convenience function to replicate the simulation results in Wu and Witten (2019). 
}
\usage{
tfCox_choose_lambda(dat, test_dat, ord = 0, alpha = 1, discrete = NULL, 
lam_seq = NULL, nlambda = 30, c = NULL, tol = 1e-06, niter=1000, 
stepSize=25, backtracking=0)
}
\arguments{
  \item{dat}{
A list that contains \code{time}, \code{status} and \code{X}. \code{time} is failure or censoring time, \code{status} is censoring indicator, and  \code{X} is n x p matrix and may have p > n. This is the training data that will be used for estimation for a given tuning parameter lambda. 
}
  \item{test_dat}{
Same list frame as before. This is the testing data that will be used for selecting tuning parameter based on the log likelihood fit. 
}
  \item{ord}{
The polynomial order of the trend filtering fit; a non-negative interger (\code{ord>= 3} is not recommended). For instance, \code{ord=0} will produce piewise constant fit, \code{ord=1} will produce piewise linear fit, and \code{ord=2} will produce piewise quadratic fit.
}
  \item{alpha}{
The trade-off between trend filtering penalty and group lasso penalty. It must be in [0,1]. \code{alpha=1} corresponds to the case with only trend filtering penalty to produce piecewise polynomial, and \code{alpha=0} corresponds to the case with only group lasso penalty to produce sparsity of the functions. \code{alpha} between 0 and 1 is the tradeoff between the strength of these two penalties. For p < n, we suggest using 1.
}
  \item{discrete}{
A vector of covariate/feature indice that are discrete. Discrete covariates are not penalized in the model. Default \code{NULL} means that none of the covariates are discrete thus all covariates will be penalized in the model.
}
  \item{lam_seq}{
The sequence of positive lambda values to consider. The default is \code{NULL}, which calculates \code{lambda.seq} using \code{lambda.min.ratio} and \code{n.lambda}. If \code{lambda.seq} is provided, it will override the default.  \code{lambda.seq} should be a decreasing positive sequence of values since \code{cv_tfCox} replies on warm starts to speed up the computation.
}
  \item{nlambda}{
The number of lambda values to consider. Default is 30.
}
  \item{c}{
Smallest value for lam_seq, as a fraction of the maximum lambda value, which is the smallest value such that the penalty term is zero. The default is NULL. 
}
  \item{tol}{
Convergence criterion for estimates.
}
  \item{niter}{
Maximum number of iterations.
}
  \item{stepSize}{
Iniitial step size. Default is 25.
}
    \item{backtracking}{
Whether backtracking should be used 1 (TRUE) or 0 (FALSE). Default is 0 (FALSE). 
}

}

\value{
  \item{lam_seq}{
  Lambda sequence considered.
}
  \item{loss}{
  Loss based on the testing data with the same length as \code{lambda.seq}
}
  \item{knots}{
  Number of knots from the training data with the same length as \code{lambda.seq}
  }
  \item{paramfit}{
  Mean square error between the estimated and true theta for the testing data. 
  }
  \item{best_lambda}{
  The lambda that achieves the minimum loss for testing data. 
  }
}

\author{
Jiacheng Wu
}

\references{
Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758
}

\seealso{
\code{\link{predict_best_lambda}}, \code{\link{negloglik}}
}

\examples{
#generate training and testing data
dat = sim_dat(n=100, zerof=0, scenario=1)
test_dat = sim_dat(n=100, zerof=0, scenario=1)

#choose the optimal tuning parameter
cv = tfCox_choose_lambda(dat, test_dat, ord=0, alpha=1)
plot(cv$lam_seq, cv$loss)

#optimal tuning parameter
cv$best_lambda

#predict the coefficients of testing covariates from the optimal tuning parameter
#from tfCox_choose_lambda object. 
theta_hat = predict_best_lambda(cv, test_dat$X)

#calculate the loss in the testing data based on the estimated coefficients theta
negloglik(test_dat, theta_hat)
}