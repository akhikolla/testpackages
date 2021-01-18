\name{cfc.survreg}
\alias{cfc.survreg}
\title{
Cause-specific competing-risk survival analysis, using parametric survival regression models
}
\description{
Convenient function to build cause-specific, parametric survival models using the \pkg{survival} package. This is followed by application of \code{cfc} function to produce cumulative incidence functions.
}
\usage{
cfc.survreg(formula, data, newdata = NULL, dist = "weibull"
  , control = survreg.control(), tout, Nmax = 100L
  , rel.tol = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Survival formula with a multi-state status variable. See \code{cfc.prepdata}.}
  \item{data}{Data frame containing variables listed in \code{formula}.}
  \item{newdata}{Data frame of structure similar to \code{data}, perhaps without the time and status columns, to be used for generating cumulative incidence curves.}
  \item{dist}{One of \code{survreg.distributions}. It can also be a vector, in which case elements \code{1} through \code{K} (number of causes) will be extracted and assigned to each cause-specific survival model. This allows for using different distributions for different causes.}
  \item{control}{List of \code{survreg} control parameters, according to \code{survreg.control}.}
  \item{tout}{Time points, along which to produce the cumulative incidence curves.}
  \item{Nmax}{Maximum number of subdivisions to be used in the \code{cfc} quadrature algorithm.}
  \item{rel.tol}{Threshold for relative error in \code{cfc} quadrature, used as a stoppage criterion. See \code{cfc} for details.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with the following elements:
\item{K}{Number of causes.}
\item{formulas}{List of formulas used in each of the \code{K} cause-specific survival regression models.}
\item{regs}{List of all cause-specific regression objects returned by \code{survreg}, one per cause. The \code{x} field of each regression object has been substituted by the model matrix from \code{newdata}.}
\item{tout}{Same as input.}
\item{cfc}{An object of class \code{cfc}, the output of applying \code{cfc} to the parametric survival regression models constructed using \code{survreg} from \pkg{survival} package.}
}

\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
}

\author{
Mansour T.A. Sharabiani, Alireza S. Mahani
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cfc.prepdata}}, \code{\link{cfc}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(bmt)
formul <- Surv(time, cause) ~ platelet + age + tcell
ret <- cfc.survreg(formul, bmt[1:300, ], bmt[-(1:300), ]
  , Nmax = 300, rel.tol = 1e-3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
