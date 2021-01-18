\name{DI}
\alias{DI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detection-Imputation algorithm
}
\description{
The Detection-Imputation algorithm computes cellwise robust estimates of the center and covariance matrix of a data set \code{X}. The algorithm alternates between the detection of cellwise outliers and their imputation combined with re-estimation of the center and covariance matrix. By default, it starts by calling \code{\link{checkDataSet}} to clean the data.
}
\usage{
DI(X, initEst = "DDCWcov", crit = 0.01, maxits = 10, quant = 0.99,
maxCol = 0.25, checkPars = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or a data frame.
  }
  \item{initEst}{An initial estimator for the center and covariance matrix. Should be one of \code{"DDCWcov"} or \code{"TSGS"}, where the latter refers to the function GSE::TSGS. The default option \code{"DDCWcov"} uses the proposal of Raymaekers and Rousseeuw (2020)  which is much faster for increasing dimension. 
  }
  \item{crit}{The algorithm converges when the subsequent estimates of the center and covariance matrix do not differ more than \code{crit} in squared Euclidean norm.
}
  \item{maxits}{Maximum number of DI-iterations.
  }
   \item{quant}{The cutoff used to detect cellwise outliers.
  }
  \item{maxCol}{The maximum number of cellwise outliers allowed in a column.
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
\item{fracNA}{
   Only retain columns and rows with fewer NAs than this fraction.
   Defaults to \eqn{0.15}.
}
  \item{\code{precScale} \cr
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
 \item{\code{silent}\cr
  Whether or not the function progress messages should be suppressed.
  Defaults to \code{FALSE}.
}
      }
   }
  }
  
\value{
 A list with components: \cr
  \itemize{
   \item{\code{center} \cr
    The final estimate of the center of the data.
    }
     \item{\code{cov} \cr
    The final estimate of the covariance matrix.
    }
    \item{\code{nits} \cr
    Number of DI-iterations executed to reach convergence.
    }
    \item{\code{Ximp} \cr
    The imputed data.
    }
    \item{\code{indcells} \cr
    Indices of the cells which were flagged in the analysis.
    }
    \item{\code{indNAs} \cr
    Indices of the NAs in the data. 
    }
     \item{\code{Zres} \cr
    Matrix with standardized cellwise residuals of the flagged cells. Contains zeroes in the unflagged cells.
    }
    \item{\code{Zres_denom} \cr
    Denominator of the standardized cellwise residuals.
    }
    \item{\code{cellPaths} \cr
    Matrix with the same dimensions as X, in which each row contains the path of least angle regression through the cells of that row, i.e. the order of the coordinates in the path (1=first, 2=second,...)
    }
    \item{\code{checkDataSet_out} \cr
    Output of the call to \code{\link{checkDataSet}} which is used to clean the data.
    }
    }
}
\references{
J. Raymaekers and P.J. Rousseeuw (2020). Handling cellwise outliers by sparse
regression and robust covariance. \emph{Arxiv: 1912.12446}. \href{https://arxiv.org/abs/1912.12446}{(link to open access pdf)}
}
\author{
J. Raymaekers and P.J. Rousseeuw
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cellHandler}}
}
\examples{
mu <- rep(0, 3)
Sigma <- diag(3) * 0.1 + 0.9
X <- MASS::mvrnorm(100, mu, Sigma)
DI.out <- DI(X)
DI.out$cov
# For more examples, we refer to the vignette:
vignette("DI_examples")
}