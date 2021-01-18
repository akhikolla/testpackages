\name{sl}
\alias{sl,Cpp_Clmbr-method}
\alias{sl}

\title{Significance Level for Changepoint}

\description{
  Significance level of a postulate value for the changepoint's x- or (x,y)-coordinates.
}

\usage{
\S4method{sl}{Cpp_Clmbr}( theta0,  method ="CLR", tolerance =0.001, output ="T" )
\S4method{sl}{Cpp_Clmbr}( theta0, alpha0,  method ="CLR", tolerance =0.001, output ="T" )
}

\arguments{
  \item{theta0}{postulate value for 'theta', the changepoint's x-coordinate.}
  \item{alpha0}{postulate value for 'alpha', the changepoint's y-coordinate.}
  \item{method}{"CLR", "MC" or "AF" which stand for conditional likelihood-ratio, conditional likelihood-ratio by Monte Carlo or approximate-F, details below.}
  \item{tolerance}{maximum absolute error in numerical integration for the "CLR" method or in Monte-Carlo evaluation for the "MC" method, not referenced for the "AF" method.}
  \item{output}{"T", "V" or "B" which stand for text, value or both.}
}

\details{
Knowles, Siegmund and Zhang (1991) reduced the conditional likelihood-ratio significance test to a probability expression based on a generic random variable.

The default method "CLR" evaluates this probability using a geometric-expectation formula that Knowles et al. also derived. This formula slightly over-estimates, but the error is negligible for significance levels below 0.20.

Method "MC" evaluates that probability expression directly by Monte Carlo simulation, which avoids the over-estimate of the "CLR" method.

Method "AF" estimates the distribution of the likelihood-ratio statistic by the related F-distribution (or chi-squared if variance is known) which would be exact for a linear model.  This method is not exact, but it is common for non-linear regression.
}

\note{
The 'tolerance' error-limit does not include the slight over-estimate that is inherent in the "CLR" method, nor the approximation inherent in the "AF" method.
}

\value{
'sl' prints-out the result but does not return a value if 'output' is "T".  'sl' returns a numeric value if 'output' is "V" or "B".
}

\examples{
#  Data for Patient B from Smith and Cook (1980)
y <- c(37.3, 47.1, 51.5, 67.6, 75.9, 73.3, 69.4, 61.5, 31.8, 19.4)
x <- 1:10
sc <- lm.br( y ~ x )

sc $ sl( 6.1 )
\donttest{sc $ sl( 6.1, 'mc' )}
\donttest{sc $ sl( 6.1, 'mc', 0.00001 )}
sc $ sl( 6.1, 88.2, 'clr' )
sc $ sl( 6.1, 88.2, 'af' )
tmp <- sc $ sl( 6.1, 88.2, 'mc', 0.001, "B" )
tmp
}

