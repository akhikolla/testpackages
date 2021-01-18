\name{predict.s2netR}
\alias{predict.s2netR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
S3 Methods for \code{\link{s2netR}} objects.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generic predict method. Wrapper for the C++ class method \code{s2net$predict}.
}
\usage{
\method{predict}{s2netR}(object, newX, type = "default", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
A \code{\link{s2netR}} object
}
  \item{newX}{
%%     ~~Describe \code{newX} here~~
A matrix with the data to make predictions. It should be in the same scale as the original data. See \code{\link{s2Data}} to see how to format the data. 
}
  \item{type}{
%%     ~~Describe \code{type} here~~
Type of predictions. One of \code{"default"} (figure it out from the train data), \code{"response", "probs", "class"}.
}
  \item{...}{
%%     ~~Describe \code{type} here~~
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A column matrix with predictions. 
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{s2netR}}, \code{\link{s2net}}
}
\examples{
data("auto_mpg")
train = s2Data(xL = auto_mpg$P1$xL, yL = auto_mpg$P1$yL,  xU = auto_mpg$P1$xU)

model = s2netR(train, 
                s2Params(lambda1 = 0.1, 
                           lambda2 = 0,
                           gamma1 = 0.1,
                           gamma2 = 100,
                           gamma3 = 0.1),
                loss = "linear",
                frame = "ExtJT",
                proj = "auto",
                fista = s2Fista(5000, 1e-7, 1, 0.8))

valid = s2Data(auto_mpg$P1$xU, auto_mpg$P1$yU, preprocess = train)
ypred = predict(model, valid$xL)
\dontrun{
if(require(ggplot2)){
  ggplot() + 
    aes(x = ypred, y = valid$yL) + geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 2)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}% use one of  RShowDoc("KEYWORDS")
