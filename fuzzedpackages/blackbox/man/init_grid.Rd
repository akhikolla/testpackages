\name{init_grid}
\alias{init_grid}
\title{
Define starting points in parameter space.
}
\description{
This function samples the space of estimated parameters. Also handles other fixed arguments that need to be passed to the function simulating the summary statistics (sample size is likely to be one such argument). The current sampling strategy is crude but achieves three desirable effects: It tries to sample the space uniformly, avoiding large gaps; it is not exactly a regular grid; and it includes replicates of some parameter points, required for good smoothing of the likelihood surface.
}
\usage{
init_grid(lower=c(par=0), upper=c(par=1), steps=NULL, 
          nUnique=NULL, nRepl=min(10L,nUnique),
          jitterFac=0.5
         )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower}{
  A vector of lower bounds for the parameters, as well as fixed arguments to be passed to the function simulating the summary statistics. Elements must be named. 
  }
  \item{upper}{
  A vector of upper bounds for the parameters, as well as fixed parameters. Elements must be named and match those of \code{lower}.  
  }
  \item{steps}{
    Number of steps of the grid, in each dimension of estimated parameters. If NULL, a default value is defined from the other arguments. If a single value is given, it is applied to all dimensions. Otherwise, this must have the same length as \code{lower} and \code{upper} and named in the same way as the variable parameters in these arguments.    
  }
  \item{nUnique}{
  Number of distinct values of parameter vectors in output. Default is an heuristic guess for good start from not too many points, computed as  \code{floor(50^((v/3)^(1/3)))} where \code{v} is the number of variable parameters.  
  }
  \item{nRepl}{
    Number of replicates of distinct values of parameter vectors in output.
  }
  \item{jitterFac}{
    Controls the amount of jitter of the points around regular grid nodes. The default value 0.5 means that a mode can move by up to half a grid step (independently in each dimension), so that two adjacent nodes moved toward each other can (almost) meet each other. 
  }
}
\value{
A data frame. Each row defines a list of arguments of vector of the function simulating the summary statistics.
}

\examples{
set.seed(123)
init_grid()
init_grid(lower=c(mu=2.8,s2=0.5,sample.size=20),
          upper=c(mu=5.2,s2=4.5,sample.size=20),
          steps=c(mu=7,s2=9),nUnique=63)
}
