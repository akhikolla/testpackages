\name{predrgasp-class}
\alias{predrgasp-class}
\docType{class}
\title{ Predictive robust GaSP class }
\description{ S4 class for the prediction of a Robust GaSP}
\section{Objects from the Class}{
Objects of this class are created and initialized with the function \code{\link{predict.rgasp}} that computes the prediction on Robust GaSP models after the Robust GaSP model has been constructed.}

\section{Slots}{
	 \describe{
	  \item{\code{call}:}{ \code{call} to \code{predict.rgasp} function where the returned object has been created.}
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

\seealso{ \code{\link{predict.rgasp}} for more details about how to make predictions based on a \code{rgasp} object.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line