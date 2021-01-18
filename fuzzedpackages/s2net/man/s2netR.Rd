\name{s2netR}
\alias{s2netR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Trains a generalized extended linear joint trained model using semi-supervised data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is a wrapper for the class \code{\link{s2net}}. It creates the C++ object and fits the model using input \code{data}.
}
\usage{
s2netR(data, params, loss = "default", frame = "ExtJT", proj = "auto", 
        fista = NULL, S3 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
A \code{\link{s2Data}} object with the (training) data.
}
  \item{params}{
%%     ~~Describe \code{params} here~~
A \code{\link{s2Params}} object with the model hyper-parameters.
}
  \item{loss}{
%%     ~~Describe \code{loss} here~~
Loss function. One of \code{"default"} (figure it out from the data), \code{"linear"} or \code{"logit"}.
}
  \item{frame}{
%%     ~~Describe \code{frame} here~~
The semi-supervised frame: \code{"ExtJT"} (the extended linear joint trained model), \code{"JT"} (the linear joint trained model from Ryan and Culp. 2015)
}
  \item{proj}{
%%     ~~Describe \code{proj} here~~
Should the unlabeled data be shifted to remove the model's effect? One of \code{"no", "yes", "auto"} (option \code{auto} shifts the unlabeled data if the angle betwen beta and the center of the data is important)
}
  \item{fista}{
%%     ~~Describe \code{fista} here~~
Fista setup parameters. An object of class \code{\link{s2Fista}}.
}
  \item{S3}{
%%     ~~Describe \code{fista} here~~
Boolean: should the method return an S3 object (default) or a C++ object?
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
Returns an object of S3 class \code{s2netR} or a C++ object of class \code{\link{s2net}}
}
\references{
%% ~put references to the literature/web site here ~
Ryan, K. J., & Culp, M. V. (2015). On semi-supervised linear regression in covariate shift problems. \emph{The Journal of Machine Learning Research, 16}(1), 3183-3217.
}
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
\code{\link{s2net}}
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
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
