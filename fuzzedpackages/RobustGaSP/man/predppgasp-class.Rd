\name{predppgasp-class}
\alias{predppgasp-class}
\docType{class}
\title{ Predicted PP GaSP class }
\description{ S4 class for the prediction of a PP GaSP model}
\section{Objects from the Class}{
Objects of this class are created and initialized with the function \code{\link{predict.ppgasp}} that computes the prediction on the PP GaSP model after the PP GaSP model has been constructed.}

\section{Slots}{
	 \describe{
	  \item{\code{call}:}{ \code{call} to \code{predict.ppgasp} function where the returned object has been created.}
    \item{\code{mean}:}{ predictive mean for the testing inputs.}
    \item{\code{lower95}:}{lower bound of the 95\% posterior credible interval.}
    \item{\code{upper95}:}{upper bound of the 95\% posterior credible interval.}
    \item{\code{sd}:}{standard deviation of each \code{testing_input}.}
  }
}


%\section{Methods}{
%  \describe{
%    \item{summary}{A summary of the object created. }
%    \item{show}{Prints the main slots of the object. }
%    \item{plot}{See \code{\link[=plot,SAVE-method]{plot}}. }
%    \item{predict}{See \code{\link[=predict.rgasp]{predict}}.}
%	 }
%}

\author{
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{predict.ppgasp}} for more details about how to make predictions based on a \code{ppgasp} object.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line