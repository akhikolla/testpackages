\name{simulate_extra}
\alias{simulate_extra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate extrapolated data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Simulated data scenarios described in the paper from Ryan and Culp (2015). 

\if{html}{\figure{sim_extra.jpg}{options: width=300 align='center'}}
}
\usage{
simulate_extra(n_source = 100, n_target = 100, p = 1000, shift = 10, 
               scenario = "same", response = "linear", sigma2 = 2.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_source}{
Number of source samples (labeled)
}
  \item{n_target}{
Number of target samples (unlabeled)
}
  \item{p}{
%%     ~~Describe \code{p} here~~
Number of variables ( \code{p > 10})
}
  \item{shift}{
%%     ~~Describe \code{shift} here~~
The shift applied to the first 10 columns of xU.
}
  \item{scenario}{
%%     ~~Describe \code{scenario} here~~
Simulation scenario. One of \code{"same"} (same distribution), \code{"lucky"} (extrapolation with lucky \eqn{\beta}), \code{"unlucky"} (extrapolation with unlucky \eqn{\beta})
}
  \item{response}{
%%     ~~Describe \code{response} here~~
Type of response: \code{"linear"} or \code{"logit"}
}
  \item{sigma2}{
%%     ~~Describe \code{sigma2} here~~
The variance of the error term, linear response case.
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
A list, with 
\describe{
  \item{xL}{ data frame with the labeled (source) data}
  \item{yL}{ labels associated with \code{xL}}
  \item{xU}{ data frame with the unlabeled (target) data}
  \item{yU}{ labels associated with \code{xU} (for validation/testing)}
}
}
\references{
%% ~put references to the literature/web site here ~
Ryan, K. J., & Culp, M. V. (2015). On semi-supervised linear regression in covariate shift problems. \emph{The Journal of Machine Learning Research, 16}(1), 3183-3217.
}
%\author{
%%  ~~who you are~~
%Juan C. Laria
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{simulate_groups}}
}
\examples{
set.seed(0)
data = simulate_extra()

train = s2Data(data$xL, data$yL, data$xU)
valid = s2Data(data$xU, data$yU, preprocess = train)

model = s2netR(train, s2Params(0.1))
ypred = predict(model, valid$xL)
plot(ypred, valid$yL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
