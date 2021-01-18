\name{buildFONKgpointls}
\alias{buildFONKgpointls}
\title{
  Prepare data for smoothing
}
\description{
  From a data frame, builds another data frame. The input data frame must contain values of the canonical parameters of the model and the variables required to construct the smoothed response. Which of the (output) parameters are variable is also determined for later use.       
}
\usage{
buildFONKgpointls(pointls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pointls}{
  A data frame obtained as return value from \code{\link{buildPointls}}
}
}
\details{
 With controls set by the Migraine software, this can operate transformations of parameter space as well as transformations in logarithmic scale (see \code{\link{islogscale}}). The output frame will then contain values of transformed parameters. 
 }
\value{
A data frame.
}
