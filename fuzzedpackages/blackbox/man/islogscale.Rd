\name{islogscale}
\alias{islogscale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for parameter log scale
}
\description{
This tests whether a log scale is used for a parameter.
}
\usage{
islogscale(string, scale = blackbox.getOption("FONKgScale"),
           extraScale = blackbox.getOption("extraScale"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{string}{
  Name of the parameter tested
}
  \item{scale}{
  A vector of scales for parameters of the smoothed object (i.e. parameters in  \code{blackbox.getOption("FONKgNames")}, see Details in \code{\link{blackbox.options}}).
}
  \item{extraScale}{
  A vector of scales for additional transformed parameters not in \code{blackbox.getOption("FONKgNames")}.
  }
}
\value{
A boolean.
}
