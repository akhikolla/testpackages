\name{BBPMM}
\alias{BBPMM}
\alias{BBPMMcol}

\title{(Multiple) Imputation through Bayesian Bootstrap Predictive Mean
  Matching (BBPMM)}

\description{\sQuote{BBPMM} performs single and multiple imputation
  (MI) of mixed-scale variables using a chained equations approach and
  (Bayesian Bootstrap) Predictive Mean Matching. 
}

\usage{BBPMM(Data, M=10, nIter=10, outfile=NULL, ignore=NULL,
vartype=NULL, stepmod="stepAIC", maxit.multi=3, maxit.glm=25,
maxPerc = 0.98, verbose=TRUE, setSeed, chainDiagnostics=TRUE, ...)}

\arguments{
  \item{Data}{A partially incomplete data frame or matrix.}
  \item{M}{Number of multiple imputations. If \env{M}=\kbd{1}, no Bayesian
    Bootstrap step is carried out. \var{Default}=\kbd{10}.}
  \item{nIter}{Number of iterations of the chained equations algorithm
  before the data set is stored as an 'imputed data set'. If set to "autolin", the numbers of iterations will be selected using a data monotonicity index (based on \code{\link[BaBooN]{dmi}}). \var{Default}=\kbd{10}.}
  \item{outfile}{A character string that specifies the path and
    file name for the imputed data sets. If \env{outfile}=\kbd{NULL} (default), no data set is
    stored}
  \item{ignore}{A character or numerical vector that specifies either
    column positions or variable names that are to be excluded from the
    imputation model and process, e.g. an ID variable. If
    \env{ignore}=\kbd{NULL} (default), all variables in \env{Data} are used
    in the imputation model.}
  \item{vartype}{A character vector that flags the class of each
    variable in \env{Data} (without the variables defined by the
    \env{ignore} argument), with either 'M' for
    metric-scale or 'C' for categorical. The default (NULL) takes over
    the classes of \env{Data}. Overruling these classes can sometimes
    make sense: e.g., an ordinal-scale variable is originally classified as \sQuote{factor},
    but treating it as metric-scale variable within the imputation
    process might still be a better choice (considering the robust properties
    of predictive mean matching to model misspecification).}
 \item{stepmod}{Performs variable selection for each imputation model based on the 
    either on Schwarz (Bayes) Information criterion (backward). \var{Default}=\kbd{"stepAIC"}.}
  \item{maxit.multi}{Imported argument from the \pkg{nnet} package that
    specifies the maximum number of iterations for the multinomial logit
    model estimation. \var{Default}=\kbd{3}.}
  \item{maxit.glm}{Argument for specification of the maximum number of iterations 
  for the binomial logit model estimation (i.e., \link[stats]{glm}). \var{Default}=\kbd{25}.}
  \item{maxPerc}{The maximum percentage the mode category of a variable is allowed to have in order to 
    try \sQuote{regular} imputation. If a variable is approximately Dirac distributed, i.e. if 
    it has (almost) no variance, imputation is carried out by simple hot deck imputation. 
    \var{Default} = \kbd{0.98}.}
  \item{verbose}{The algorithm prints information on imputation and
    iteration numbers. \var{Default}=\kbd{TRUE}.}
  \item{setSeed}{Optional argument to fix the pseudo-random number generator in order to
    allow for reproducible results.}
  \item{chainDiagnostics}{
    Argument specifying if Monte Carlo chains for further diagnostics should be returned as well. \var{Default}=\kbd{TRUE}.}
  \item{...}{Further arguments passed to or from other functions.}
}

\value{
  \item{call}{The call of \command{BBPMM}.}
  \item{mis.num}{Vector containing the numbers of missing values per column.}
  \item{modelselection}{Chosen model selection method for the function call.}
  \item{seed}{Chosen seed value for the function call.}
  \item{impdata}{The imputed data set, if \env{M}=\kbd{1}, or a list containing \var{M}
    imputed data sets.}
  \item{misOverview}{The percentage of missing values per incomplete variable.}
  \item{indMatrix}{A matrix with the same dimensions as \env{Data}
    minus \env{ignore} containing flags for missing values.}
  \item{M}{Number of (multiple) imputations.}
  \item{nIter}{Number of iterations between two imputations.}
  \item{Chains}{List containing the the Gibbs sampler sequences for every variable of every imputation for every iteration.}
  \item{FirstSeed}{First \code{\link[=Random]{.Random.Seed}} before imputation starts.}
  \item{LastSeed}{Last \code{\link[=Random]{.Random.Seed}} after function is done.}
  \item{ignoredvariables}{\kbd{TRUE} / \kbd{FALSE} indicator whether variables were ignored during imputation.}
}

\details{\code{BBPMM} is based on a chained equations approach
  that is using a Bayesian Bootstrap approach and Predictive Mean
  Matching (\acronym{PMM}) variants for metric-scale, binary, and multi-categorical
  variables to generate multiple imputations. In order to emulate a
  monotone missing-data pattern as well as possible, variables are sorted
  by rate of missingness (in ascending order). If no complete variables
  exist, the least incomplete variable is imputed via hot-deck. The
  starting solution then builds the imputation model using the observed values of
  a particular y variable, and the corresponding observed or already
  imputed values of the \var{x} variables (i.e., all variables with fewer
  missing values than \var{y}).
  Due to the \acronym{PMM} element in the algorithm,
  auto-correlation of subsequent iterations is virtually zero. Therefore, a
  burn-in period is not required, and there is no need to administer
  \sQuote{high} values (> 20) to \option{nIter} either.
  
  If \env{M}=\kbd{1}, no Bayesian Bootstrap step is carried
  out for the chained equations. Note that in this case the algorithm is still unlikely to
  converge to a stable solution, because of the Predictive Mean Matching
  step.
}

\author{Florian Meinfelder, Thorsten Schnapp [ctb]}


\seealso{\code{\link{BBPMM.row}}, \code{\link{dmi}}}

\references{

  Koller-Meinfelder, F. (2009) \emph{Analysis of Incomplete Survey Data
    -- Multiple Imputation Via Bayesian Bootstrap Predictive Mean
    Matching}, doctoral thesis.

  Little, R.J.A. (1988) Missing-Data Adjustments in Large
    Surveys, \emph{Journal of Business and Economic Statistics}, Vol. \bold{6},
    No. 3, pp. 287-296.

  Raghunathan T.E. and Lepkowski, J.M. and Van Hoewyk, J. and Solenberger, P (2001) A multivariate technique
    for multiply imputing missing values using a sequence of regression models. \emph{Survey Methodology}, Vol. \bold{27}, pp. 85--95.

  Rubin DB (1981) The Bayesian Bootstrap. \emph{The Annals of Statistics}, Vol. \bold{9}, pp. 130--134.

	 Rubin, D.B. (1987) \emph{Multiple Imputation for Non-Response in Surveys}. New York: John Wiley & Sons, Inc.

  Van Buuren, S. and Brand, J.P.L. and Groothuis-Oudshoorn, C.G.M. and Rubin, D.B. (2006) Fully conditional specification
  in multivariate imputation. \emph{Journal of Statistical Computation and Simulation}, Vol. \bold{76}, No. 12, pp. 1049--1064.

  Van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation by Chained Equations in R. \emph{Journal of
    Statistical Software}, Vol. \bold{45}, No. 3, pp. 1--67. URL http://www.jstatsoft.org/v45/i03/.

  Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth Edition. New York: Springer.

}

\keyword{robust}
\keyword{datagen}

\examples{

### sample data set with non-normal variables
set.seed(1000)
n <- 50
x1 <- round(runif(n,0.5,3.5))
x2 <- as.factor(c(rep(1,10),rep(2,25),rep(3,15)))
x3 <- round(rnorm(n,0,3))
y1 <- round(x1-0.25*(x2==2)+0.5*x3+rnorm(n,0,1))
y1 <- ifelse(y1<1,1,y1)
y1 <- as.factor(ifelse(y1>4,5,y1))
y2 <- x1+rnorm(n,0,0.5)
y3 <- round(x3+rnorm(n,0,2))
data1 <- as.data.frame(cbind(x1,x2,x3,y1,y2,y3))
misrow1 <- sample(n,20)
misrow2 <- sample(n,15)
misrow3 <- sample(n,10)
is.na(data1[misrow1, 4]) <- TRUE
is.na(data1[misrow2, 5]) <- TRUE
is.na(data1[misrow2, 6]) <- TRUE

### imputation
imputed.data <- BBPMM(data1, nIter=5, M=5)
}
