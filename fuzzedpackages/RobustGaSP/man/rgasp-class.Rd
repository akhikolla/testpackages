\name{rgasp-class}
\docType{class}
\alias{rgasp-class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust GaSP class }
\description{ S4 class for Robust GaSP if the range and noise-variance ratio parameters are given and/or have been estimated.}

\section{Objects from the Class}{
Objects of this class are created and initialized with the function \code{\link{rgasp}} that computes the calculations needed for setting up the analysis.}

\section{Slots}{
	 \describe{
    \item{\code{p}:}{Object of class \code{integer}. The dimensions of the inputs.}
    \item{\code{num_obs}:}{Object of class \code{integer}. The number of observations.}
    \item{\code{input}:}{Object of class \code{matrix} with dimension n x p. The design of experiments.}
    \item{\code{output}:}{Object of class \code{matrix} with dimension n x 1. The Observations or output vector.}
    \item{\code{X}:}{Object of class \code{matrix} of with dimension n x q. The mean basis function, i.e. the trend function.}
        \item{\code{zero_mean}:}{A \code{character} to specify whether the mean is zero or not. "Yes" means it has zero mean and "No"" means the mean is not zero. }
    \item{\code{q}:}{Object of class \code{integer}. The number of mean basis.}
    \item{\code{LB}:}{Object of class \code{vector} with dimension p x 1. The lower bound for inverse range parameters beta.}
    \item{\code{beta_initial}:}{Object of class \code{vector} with the initial values of inverse range parameters p x 1.}
    \item{\code{beta_hat}:}{Object of class \code{vector} with dimension p x 1. The inverse-range parameters.}
    \item{\code{log_post}:}{Object of class \code{numeric} with the logarithm of marginal posterior.}
    \item{\code{R0}:}{Object of class \code{list} of matrices where the j-th matrix is an absolute difference matrix of the j-th input vector.}
    \item{\code{theta_hat}:}{Object of class \code{vector} with dimension q x 1. The the mean (trend) parameter.}
    \item{\code{L}:}{Object of class \code{matrix} with dimension n x n. The Cholesky decomposition of the correlation matrix \code{R}, i.e. \deqn{L\%*\%t(L)=R}{L\%*\%t(L)=R}}
    \item{\code{sigma2_hat}:}{Object of the class \code{numeric}. The estimated variance parameter.}
    \item{\code{LX}:}{Object of the class \code{matrix} with dimension q x q. The Cholesky decomposition of the correlation matrix \deqn{t(X)\%*\%R^{-1}\%*\%X}{t(X)\%*\%R^{-1}\%*\%X}}
    \item{\code{CL}:}{Object of the class \code{vector} used for the lower bound and the prior.}
    \item{\code{nugget}:}{A \code{numeric} object used for the noise-variance ratio parameter.}
    \item{\code{nugget.est}:}{A \code{logical} object of whether the nugget is estimated (T) or fixed (F).}
    \item{\code{kernel_type}:}{A \code{vector} of \code{character} to specify the type of kernel to use.}
    \item{\code{alpha}:}{Object of class \code{vector} with dimension p x 1 for the roughness parameters in the kernel.}
        \item{\code{method}:}{Object of class \code{character} to specify the method of parameter estimation. There are three values: \code{post_mode}, \code{mle} and \code{mmle}.}
        \item{\code{isotropic}:}{Object of class \code{logical} to specify whether the kernel is isotropic. }
    \item{\code{call}:}{The \code{call} to \code{rgasp} function to create the object.}
    }
  }


\section{Methods}{
  \describe{
%    \item{summary}{A summary of the object created. }
    \item{show}{Prints the main slots of the object. }
%    \item{plot}{See \code{\link[=plot,SAVE-method]{plot}}. }
    \item{predict}{See \code{\link[=predict.rgasp]{predict}}.}
	 }
}

\author{
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

\note{
%%  ~~further notes~~
The response \code{output} must have one dimension.
The number of observations in \code{input} must be equal to the number of experiments \code{output}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{RobustGaSP}} for more details about how to create a \code{RobustGaSP} object.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
