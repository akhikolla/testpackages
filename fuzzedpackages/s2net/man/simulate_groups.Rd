\name{simulate_groups}
\alias{simulate_groups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate data (two groups design)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Simulated data scenario described in paper [citation here].

\if{html}{\figure{sim_fr.jpg}{options: width=600 align='center'}}
}
\usage{
simulate_groups(n_source = 100, n_target = 100, p = 200, response = "linear")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_source}{
%%     ~~Describe \code{n_source} here~~
Number of labeled observations
}
  \item{n_target}{
%%     ~~Describe \code{n_target} here~~
Number of unlabeled (target) observations
}
  \item{p}{
%%     ~~Describe \code{p} here~~
Number of variables
}
  \item{response}{
%%     ~~Describe \code{response} here~~
Type of response: \code{"linear"} or \code{"logit"}
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
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
%%  ~~who you are~~
Juan C. Laria
}
%\note{
%%  ~~further notes~~
%}
%
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{simulate_extra}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
