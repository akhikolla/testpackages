\name{warpNorm}
\alias{warpNorm}
\alias{warpPower}
%\docType{data}
\title{
Warpings for ordinal inputs
}
\description{
Given warpings for ordinal inputs.
}
\usage{
warpNorm
warpPower
}
\format{
  The format is a list of 6:

 $ fun : the warping function. The second argument is the vector of
 parameters. The function returns a numeric vector with an attribute
 \code{"gradient"} giving the derivative w.r.t. the parameters.
 
 $ parNames  : names of warping parameters (character vector). 
 
 $ parDefault: default values of warping parameters (numeric vector).
 
 $ parLower  : lower bounds of warping parameters (numeric vector).
 
 $ parUpper  : upper bounds of warping parameters (numeric vector). 
 
 $ hasGrad : a boolean equal to \code{TRUE} if \code{gradient} is
 supplied as an attribute of \code{fun}.

} 
\details{
  
  See \code{\link{covOrd}} for the definition of a warping in this
  context.  At this stage, two warpings corresponding to cumulative
  density functions (cdf) are implemented:
  
  \itemize{
    
    \item Normal distribution, truncated to \eqn{[0,1]}: \deqn{F(x) = [N(x) -
      N(0)] / [N(1) - N(0)]} where \eqn{N(x) = \Phi([x - \mu] /
      \sigma)}{N(x) = Phi([x - mu] / sigma)} is the cdf of the normal
    distribution with mean \eqn{\mu}{mu} and standard deviation
    \eqn{\sigma}{sigma}.
    
    \item Power distribution on \eqn{[0, 1]}: \eqn{F(x) = x^{pow}}{F(x) = x^pow}.

  }
  

}
