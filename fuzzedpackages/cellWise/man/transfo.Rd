\name{transfo}
\alias{transfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robustly fit the Box-Cox or Yeo-Johnson transformation 
}
\description{
This function uses reweighted maximum likelihood to robustly fit the 
Box-Cox or Yeo-Johnson transformation to each variable in a dataset. 
Note that this function first calls \code{\link{checkDataSet}} to ensure that the variables to be transformed are not too discrete.
}
\usage{
transfo(X, type = "YJ", robust = TRUE, lambdarange = NULL,
        prestandardize = TRUE, prescaleBC = F, scalefac = 1,
        quant = 0.99, nbsteps = 2, checkPars = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ A data matrix of dimensions n x d.
  Its columns are the variables to be transformed.
}
\item{type}{ The type of transformation to be fit. Should be one of: \itemize{
\item{\code{"BC"}}: Box-Cox power transformation. Only works
 for strictly positive variables. If this type
 is given but a variable is not strictly
  positive, the function stops with a
   message about that variable.
\item{\code{"YJ"}} Yeo-Johnson power transformation. The data
 may have positive as well as negative values.
\item{\code{"bestObj"}} for strictly positive variables both BC
                         and YJ are run, and the solution with
   lowest objective is kept. On the other
  variables YJ is run.
}

}
\item{robust}{ if \code{TRUE} the Reweighted Maximum 
                Likelihood method is used, which first
   computes a robust initial estimate of the
 transformation parameter lambda. If \code{FALSE} the classical ML method is used.
 }
 \item{lambdarange}{ range of lambda values that will be optimized over. 
      If \code{NULL}, the range goes from -4 to 6.
      }
 \item{prestandardize}{ whether to standardize the variables \strong{before} the
 power transformation.For BC the variable is divided by its median.
For YJ and \code{robust = TRUE} this subtracts its median and divides by its mad (median absolute deviation). For YJ and \code{robust = F} this subtracts the mean and divides by the standard deviation.
 }
 \item{prescaleBC}{ for BC only. This standardizes the logarithm of the
   original variable by subtracting its median and
  dividing by its mad, after which the exponential function turns
  the result into a positive variable again.
  }
  \item{scalefac}{when YJ is fit and \code{prestandardize = TRUE}, the 
  standardized data is multiplied by \code{scalefac}. 
  When BC is fit and \code{prescaleBC = TRUE} the same happens 
  to the standardized log of the original variable.
  }
 \item{quant}{quantile for determining the weights in the 
   reweighting step (ignored when \code{robust=FALSE}).
   }
 \item{nbsteps}{number of reweighting steps (ignored when 
   \code{robust=FALSE}).
   }
 \item{checkPars}{Optional list of parameters used in the call to
 \code{\link{checkDataSet}}. The options are:
  \itemize{
  \item{ \code{coreOnly} \cr
      If \code{TRUE}, skip the execution of checkDataset. Defaults to \code{FALSE}}
       \item{\code{numDiscrete}\cr
  A column that takes on numDiscrete or fewer values
  will be considered discrete and not retained in the cleaned data.
  Defaults to \eqn{5}.
}
  \item{\code{precScale} \cr
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
 \item{\code{silent}\cr
  Whether or not the function progress messages should be printed.
  Defaults to \code{FALSE}.
}
      }
   }
}
\value{
A list with components: \cr
  \itemize{
    \item{\code{lambdahats} \cr
    the estimated transformation parameter for each column of \code{X}.
    }
    \item{\code{Xt} \cr
    A matrix in which each column is the transformed version of the
   corresponding column of \code{X}.
    }
      \item{\code{muhat} \cr
    The estimated location of each column of \code{Xt}.
    }
      \item{\code{sigmahat} \cr
    The estimated scale of each column of \code{Xt}.
    }
     
     \item{\code{Zt} \cr
    \code{Xt} poststandardized by the centers in \code{muhat} and the
    scales in \code{sigmahat}. Is always provided.
    }
     \item{\code{weights} \cr
    The final weights from the reweighting.
    }
      \item{\code{ttypes} \cr
    The type of transform used in each column.
    }
     \item{\code{objective} \cr
     Value of the (reweighted) maximum likelihood objective function.
    }
              
    }
}
\references{
J. Raymaekers and P.J. Rousseeuw (2020). Transforming variables to 
central normality. \emph{Arxiv: 2005.07946}. \href{https://arxiv.org/abs/2005.07946}{(link to open access pdf)}
}
\author{
J. Raymaekers and P.J. Rousseeuw
}

\examples{

# find Box-Cox transformation parameter for lognormal data:
set.seed(123)
x <- exp(rnorm(1000))
transfo.out <- transfo(x, type = "BC")
# estimated parameter:
transfo.out$lambdahat
# value of the objective function:
transfo.out$objective
# the transformed variable:
transfo.out$Xt
# the poststandardized transformed variable:
transfo.out$Zt
# the type of transformation used:
transfo.out$ttypes
# qqplot of the poststandardized transformed variable:
qqnorm(transfo.out$Zt); abline(0,1)

# For more examples, we refer to the vignette:
vignette("transfo_examples")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
