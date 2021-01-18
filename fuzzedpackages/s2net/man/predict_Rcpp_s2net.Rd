\name{predict_Rcpp_s2net}
\alias{predict_Rcpp_s2net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Predict method for \code{s2net} C++ class.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function provides an interface in \R for the method \code{predict} in C++ class \code{s2net}.
}
\usage{
%% \method{predict}{Rcpp_s2net}(object, newX, type = "default", ... )
predict_Rcpp_s2net(object, newX, type = "default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
An object of class \code{\link{Rcpp_s2net}}.
}
  \item{newX}{
%%     ~~Describe \code{newX} here~~
Data to make predictions. Could be a \code{\link{s2Data}} object (field \code{xL} is used) or a \code{matrix} (in the same space as the original data where the model was fitted).  
}
  \item{type}{
%%     ~~Describe \code{type} here~~
Type of predictions. One of \code{"default"}: let the method figure it out; \code{"response"}: the linear predictor; \code{"probs"}: fitted probabilities; \code{class}: fitted class.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This method is included as a high-level wrapper of \code{object$predict()}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a column \code{matrix} with the same number of rows/observations as \code{newX}.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
%%  ~~who you are~~
Juan C. Laria
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{Rcpp_s2net}}
}
%%\examples{

%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{}% use one of  RShowDoc("KEYWORDS")
