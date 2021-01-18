\name{log_approx_ref_prior_deriv}
\alias{log_approx_ref_prior_deriv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Derivative of the jointly robust prior
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function computes the derivative of the approximate reference prior with regard to inverse range parameter and the nugget-noise ratio parameter (if not fixed). When the nugget is fixed, it only compute the derivative with regard to the inverse range parameter; otherwise it produces derivative with regard to inverse range parameter and noise ratio parameter.
}
\usage{
log_approx_ref_prior_deriv(param, nugget, nugget_est, CL, a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
A vector of natural logarithm of inverse-range parameters and natural logarithm of the nugget-variance ratio parameter.
}
  \item{nugget}{
%%     ~~Describe \code{nugget} here~~
The nugget-variance ratio parameter if this parameter is fixed.
}
  \item{nugget_est}{
%%     ~~Describe \code{nugget.est} here~~
Boolean value of whether the nugget is estimated or fixed.
}
  \item{CL}{
%%     ~~Describe \code{CL} here~~
Prior parameter in the jointly robust prior.
}
  \item{a}{
%%     ~~Describe \code{a} here~~
Prior parameter in the jointly robust prior.
}
  \item{b}{
%%     ~~Describe \code{b} here~~
Prior parameter in the jointly robust prior.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The numerical value of the derivative of the jointly robust prior with regard to \code{beta} (the inverse-range parameters) and \code{nugget} (the nugget-variance ratio parameter). When the nugget is fixed, the derivative is on inverse-range parameters.
}
\references{
M. Gu. (2016). Robust uncertainty quantification and scalable computation for computer models with massive output. Ph.D. thesis. Duke University.
}
\author{
%%  ~~who you are~~
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{log_approx_ref_prior}},\code{\link{rgasp}} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}