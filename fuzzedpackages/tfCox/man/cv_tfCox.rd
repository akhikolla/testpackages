\name{cv_tfCox}
\alias{cv_tfCox}
\title{Fit Trend Filtering Cox model and Choose Tuning Parameter via K-Fold Cross-Validation}
\description{Fit additive trend filtering Cox model where each component function is estimated to be piecewise constant or polynomial. Tuning parameter is selected via k-fold cross-validation. }
\usage{
cv_tfCox(dat, ord=0, alpha=1, discrete=NULL, lambda.seq=NULL,
lambda.min.ratio=0.01, n.lambda=30, n.fold=5, seed=NULL, tol=1e-6,
niter=1000, stepSize=25,  backtracking=0)
}
\arguments{
  \item{dat}{
A list that contains \code{time}, \code{status} and \code{X}. \code{time} is failure or censoring time, \code{status} is censoring indicator, and  \code{X} is n x p matrix and may have p > n.
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
  \item{lambda.seq}{
The sequence of positive lambda values to consider. The default is \code{NULL}, which calculates \code{lambda.seq} using \code{lambda.min.ratio} and \code{n.lambda}. If \code{lambda.seq} is provided, it will override the default.  \code{lambda.seq} should be a decreasing positive sequence of values since \code{cv_tfCox} replies on warm starts to speed up the computation.
}
  \item{lambda.min.ratio}{
Smallest value for lambda.seq, as a fraction of the maximum lambda value, which is the smallest value such that the penalty term is zero. The default is 0.01.
  }
  \item{n.lambda}{
The number of lambda values to consider. Default is 30.
}
  \item{n.fold}{
The number of folds for cross-validation of \code{lambda}. The default is 5.
  }
  \item{seed}{
An optional number used with \code{set.seed()}.
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

\details{
Note that \code{cv_tfCox} does not cross-validate over \code{alpha}, and \code{alpha} should be provided. However, if the user would like to cross-validate over \code{alpha}, then \code{cv_tfCox} should be called multiple times for different values of \code{alpha} and the same \code{seed}. This ensures that the cross-validation folds (\code{fold}) remain the same for the different values of \code{alpha}. See the example below for details.
}

\value{
An object with S3 class "cv_tfCox".
  \item{best.lambda}{
Optional lambda value chosen by cross-dalidation.
}
  \item{lambda.seq}{
lambda sequence considered.
}
  \item{mean.cv.error}{
vector of average cross validation error with the same length as \code{lambda.seq}
}
}

\author{
Jiacheng Wu
}

\references{
Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758
}

\seealso{
\code{\link{summary.cv_tfCox}}, \code{\link{plot.cv_tfCox}}, \code{\link{tfCox}}
}

\examples{
#generate data
set.seed(123)
dat = sim_dat(n=100, zerof=0, scenario=1)

#fit piecewise constant functions
#cross-validation to choose the tuning parameter lambda with fixed alpha=1
cv = cv_tfCox(dat, ord=0, alpha=1, n.fold=2, seed=123)
plot(cv, showSE=TRUE)
}
