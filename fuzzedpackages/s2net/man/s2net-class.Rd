\name{Rcpp_s2net-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_s2net-class}
\alias{Rcpp_s2net}
\alias{_rcpp_module_boot_Rcpp_s2net_export}
\alias{predict,Rcpp_s2net-method}

\title{Class \code{s2net}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
This is the main class of this library, implemented in C++ and exposed to \R using \code{Rcpp} modules. 
It can be used in \R directly, although some generic S4 methods have been implemented to make it easier to interact in \R.
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Rcpp_s2net")}: See \code{\link{predict_Rcpp_s2net}}}
	 }
}
%%\references{
%%  ~~put references to the literature/web site here~~
%%}
\author{
%%  ~~who you are~~
Juan C. Laria
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%%}
\examples{
data("auto_mpg")
train = s2Data(xL = auto_mpg$P1$xL, yL = auto_mpg$P1$yL,  xU = auto_mpg$P1$xU)

# We create the C++ object calling the new method (constructor)
obj = new(s2net, train, 0) # 0 = regression 
obj

# We call directly the $fit method of obj, 
obj$fit(s2Params(lambda1 = 0.01, 
                   lambda2 = 0.01, 
                   gamma1 = 0.05, 
                   gamma2 = 100, 
                   gamma3 = 0.05), 1, 2)
# fitted model
obj$beta

# We can test the results using the unlabeled data
test = s2Data(xL = auto_mpg$P1$xU, yL = auto_mpg$P1$yU,  preprocess = train)
ypred = obj$predict(test$xL, 0)

\dontrun{
if(require(ggplot2)){
  ggplot() + 
    aes(x = ypred, y = test$yL) + geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 2)
}
}

}
\keyword{optimize}
\keyword{models}
\keyword{regression}

\section{Fields}{
  \describe{
    \item{\code{beta}:}{Object of class \code{matrix}. The fitted model coefficients.}
    \item{\code{intercept}:}{ The model intercept.}
  }
}
\section{Class-Based Methods}{
  \describe{
    \item{\code{initialize(data, loss)}:}{
      \describe{
      \item{\code{data} }{ \code{\link{s2Data}} object}
      \item{\code{loss} }{Loss function: 0 = linear, 1 = logit}
      }
    }
    \item{\code{setupFista(\link{s2Fista})}:}{ Configures the FISTA internal algorithm.}
    \item{\code{predict(newX, type)}:}{
      \describe{
      \item{\code{newX} }{ New data \code{matrix} to make predictions.}
      \item{\code{type} }{ 0 = default, 1 = response, 2 = probs, 3 = class}
      }
    }
    \item{\code{fit(params, frame, proj)}:}{
    \describe{
      \item{\code{params }}{\code{\link{s2Params}} object}
      \item{\code{frame} }{0 = "JT", 1 = "ExtJT"}
      \item{\code{proj} }{0 = no, 1 = yes, 2 = auto}
    }
    }
  }
}
