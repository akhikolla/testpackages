
\name{'lm.br'}
\alias{lm.br-package}
\alias{lm.br}
\alias{print.lm.br}


\docType{package}

\title{Fit a Linear Model with a Breakpoint}

\description{
Exact significance tests for a changepoint in
linear or multiple linear regression.  Confidence intervals
and confidence regions with exact coverage
probabilities for the changepoint.
}


\usage{
lm.br(formula, type ="LL", data, subset, weights, inverse =FALSE,
      var.known =FALSE, na.action, contrasts, offset, ...)
}


\arguments{
  \item{formula}{a formula expression as for regression models, of
the form \code{response ~ predictors}; see \code{\link{formula}}.}
  \item{type}{ "LL", "LT" or "TL" which stand for line-line,
line-threshold or threshold-line, defined below.}
  \item{data}{an optional data-frame that assigns values in
\code{formula}.}
  \item{subset}{expression saying which subset of the data to use.}
  \item{weights}{vector or matrix.}
  \item{inverse}{if TRUE then 'weights' specifies the inverse of the
weights vector or matrix, as for a covariance matrix.}
  \item{var.known}{is the variance known?}
  \item{na.action}{a function to filter missing data.}
  \item{contrasts}{an optional list; see 'contrasts.arg' in
\code{\link{model.matrix}}.}
  \item{offset}{a constant vector to be subtracted from the
responses vector.}
  \item{\dots}{ other arguments to \code{\link{lm.fit}} or
\code{\link{lm.wfit}}.}
}


\details{
  A broken-line model consists of two straight lines joined at a
  changepoint.  Three versions are
\preformatted{LL   y = alpha + B * min(x - theta, 0) + Bp * max(x - theta, 0) + e

LT   y = alpha + B * min(x - theta, 0) + e

TL   y = alpha + Bp * max(x - theta, 0) + e}
  where  e ~ Normal( 0, var * inv(weights) ).  The LT and TL versions
  omit 'alpha' if the formula is without intercept, such as 'y~x+0'.  
  Parameters 'theta', 'alpha', 'B', 'Bp', 'var' are unknown, but 
  'weights' is known.
  
  The same models apply for a multiple-regression formula such as 'y ~ x1 +
  x2 + ... + xn'  where 'alpha' becomes the coefficient of the
  "1"-vector and 'theta' the changepoint for the coefficient of the 
  first predictor term, 'x1'.

  The test for the presence of a changepoint is by
  a postulate value outside the range of 'x'-values.  Thus, in the 
  LL model 'sl( min(x1) - 1 )'  would give the exact significance 
  level of the null hypothesis "single line" versus the alternate 
  hypothesis "broken line."

  Exact inferences about the changepoint  
  'theta'  or '(theta,alpha)'  are based on the distribution of its 
  likelihood-ratio statistic, conditional on sufficient statistics 
  for the other parameters.  This method is called conditional likelihood-ratio (CLR) for short.
}


\value{
'lm.br' returns a list that includes a C++ object with accessor
functions.  Functions \code{\link{sl}}, \code{\link{ci}} and \code{\link{cr}} get significance levels, confidence intervals, 
and confidence regions for the changepoint's x-coordinate or
(x,y)-coordinates.  Other functions are \code{\link{mle}} to get maximum likelihood estimates and \code{\link{sety}} to set new y-values.  
The returned object also lists 'coefficients', 'fitted.values' and 'residuals', the same as for an 'lm' output list.
}


\note{
Data can include more than one 'y' value for a repeat 'x' value.  
If variance is known, then  'var' = 1  and  'weights'  is the inverse of the variances
vector or variance-covariance matrix.
}


\references{
Knowles, M., Siegmund, D. and Zhang, H.P. (1991)  Confidence regions 
in semilinear regression,  _Biometrika_, *78*, 15-31. 

Siegmund, D. and Zhang, H.P. (1994),  Confidence regions in
broken line regression,  in  "Change-point Problems",  _IMS
Lecture Notes -- Monograph Series_, *23*,  eds. E. Carlstein, H.
Muller and D. Siegmund,  Hayward, CA: Institute of Mathematical
Statistics,  292-316.
}

\keyword{ package }

\seealso{
vignette( "\href{../doc/lm.br.pdf}{lm.br}" )
\cr demo( testscript )
}

\examples{
#  Smith & Cook (1980), "Straight Lines with a Change-point: A Bayesian
#  Analysis of some Renal Transplant Data", Appl Stat, *29*, 180-189,  
#  reciprocal of blood creatinine L/micromol  vs  day after transplant.
creatinine <- c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
day <- 1:10
sc <- lm.br( creatinine ~ day )
sc $ mle()
sc $ ci()
sc $ sl( day[1] - 1.5 )      # test for the presence of a changepoint
plot( sc$residuals )


#  A 'TL' example, data from figure 1 in Chiu et al. (2006), "Bent-cable
#  regression theory and applications", J Am Stat Assoc, *101*, 542-553,
#  log(salmon abundance) vs year.
salmon <- c( 2.50, 2.93, 2.94, 2.83, 2.43, 2.84, 3.06, 2.97, 2.94, 2.65,
  2.92, 2.71, 2.93, 2.60, 2.12, 2.08, 1.81, 2.45, 1.71, 0.55, 1.30 )
year <- 1980 : 2000
chiu <- lm.br( salmon ~ year, 'tl' )
chiu $ ci()


#  A multiple regression example, using an R dataset,
#  automobile miles-per-gallon  versus  weight and horsepower.
lm.br( mpg ~ wt + hp,  data = mtcars )


#  An example with variance known, for the Normal approximations of binomial
#  random variables using formula 2.28 of Cox and Snell (1989).
#    Ex. 3.4 of Freeman (2010) "Inference for binomial changepoint data" in
# _Advances in Data Analysis_, ed. C Skiadas, Boston: Birkhauser, 345-352.
trials <- c( 15, 82, 82, 77, 38, 81, 12, 97, 33, 75,
  85, 37, 44, 96, 76, 26, 91, 47, 41, 35 )
successes <- c( 8, 44, 47, 39, 24, 38, 3, 51, 16, 43,
  47, 27, 33, 64, 41, 18, 61, 32, 33, 24 )
log_odds <- log( (successes - 0.5)/(trials - successes - 0.5) )
variances <- (trials-1)/( successes*(trials-successes) )
group <- 1 : 20
lm.br( log_odds ~ group, 'TL', w= variances, inv= TRUE, var.known= TRUE )


#  An example that shows different confidence regions from inference by
#  conditional likelihood-ratio (CLR)  versus  approximate-F (AF).
y <- c( 1.6, 3.2, 6.3, 4.8, 4.3, 4.0, 3.5, 1.8 )
x <- 1:8
eg <- lm.br( y ~ x )
\donttest{eg$cr( output='t' )}
eg$cr( method = 'aF', output='t' )
}

