\name{predict.gp}
\alias{predict.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prediction Method for the \code{"gp"} S3 Class
}
\description{
  
  Prediction method for the \code{"gp"} S3 class.
 
}
\usage{

\S3method{predict}{gp}(object, newdata,
        type = ifelse(object$trendKnown, "SK", "UK"), 
        seCompute = TRUE, covCompute = FALSE,
        lightReturn = FALSE, biasCorrect = FALSE,
        forceInterp,        
        ...)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    An object with S3 class \code{"gp"}.

  }
  \item{newdata}{

    A data frame containing all the variables required for prediction: inputs and trend variables, if applicable.

  }
  \item{type}{

    A character string corresponding to the GP "kriging" family, to be chosen between simple kriging (\code{"SK"}), or universal kriging
    (\code{"UK"}).

  }
  \item{seCompute}{

    Optional logical. If \code{FALSE}, only the kriging mean is computed.  If \code{TRUE}, the kriging variance (actually, the corresponding standard deviation) and prediction intervals are computed too.

  }
  \item{covCompute}{

    Logical. If \code{TRUE} the covariance matrix is computed.

  }
  \item{lightReturn}{

    Optional logical. If \code{TRUE}, \code{c} and \code{cStar} are
    not returned. This should be reserved to expert users who want to
    save memory and know that they will not miss these values.

  }
  \item{biasCorrect}{

    Optional logical to correct bias in the UK variance and
    covariances. Default is \code{FALSE}. See \bold{Details}
    below.

  }
  \item{forceInterp}{
    
    Logical used to force a nugget-type prediction. If \code{TRUE},
    the noise will be interpreted as a nugget effect. \emph{This argument
    is likely to be removed in the future}.
    
  }
  \item{\dots}{

    Not used yet.

  }
}

\details{

  The estimated (UK) variance and covariances are NOT multiplied by
  \eqn{n/(n-p)} by default (\eqn{n} and \eqn{p} denoting the number of
  rows and columns of the trend matrix \eqn{\mathbf{F}}{F}). Recall that
  this correction would contribute to limit bias: it would totally
  remove it if the correlation parameters were known (which is not the
  case here). However, this correction is often ignored in the context
  of computer experiments, especially in adaptive strategies. It can be
  activated by turning \code{biasCorrect} to \code{TRUE}, when
  \code{type = "UK"}

}

\value{

  A list with the following elements.
  
  \item{mean}{

    GP mean ("kriging") predictor (including the trend) computed at
    \code{newdata}.

  }
  \item{sd}{

    GP prediction ("kriging") standard deviation computed at
    \code{newdata}. Not computed if \code{seCompute} is \code{FALSE}.

  }
  \item{sdSK}{

    Part of the above standard deviation corresponding to simple kriging
    (coincides with \code{sd} when \code{type = "SK"}). Not computed if
    \code{seCompute} is \code{FALSE}.

  }
  \item{trend}{

    The computed trend function, evaluated at \code{newdata}.

  }
  \item{cov}{

    GP prediction ("kriging") conditional covariance matrix. Not
    computed if \code{covCompute} is \code{FALSE} (default).

  }
  \item{lower95,}{ }
  \item{upper95}{

    Bounds of the 95 \% GP prediction interval computed at
    \code{newdata} (to be interpreted with special care when parameters
    are estimated, see description above). Not computed if
    \code{seCompute} is \code{FALSE}.

  }
  \item{c}{

    An auxiliary matrix \eqn{\mathbf{c}}{c}, containing all the
    covariances between the points in \code{newdata} and those in the
    initial design. Not returned if \code{lightReturn} is
    \code{TRUE}.

  }
  \item{cStar }{

    An auxiliary vector, equal to \eqn{\mathbf{L}^{-1}\mathbf{c}}{L^(-1)
    \%*\% c} where \eqn{\mathbf{L}}{L} is the Cholesky root of the
    covariance matrix \eqn{\mathbf{C}}{C} used in the estimation. Not
    returned if \code{lightReturn} is \code{TRUE}.

  }
}

\author{
  O. Roustant, D. Ginsbourger, Y. Deville
}

\seealso{\code{\link{gp}} for the creation/estimation of a model. See
  \code{\link{gls-methods}} for the signification of the auxiliary variables.}