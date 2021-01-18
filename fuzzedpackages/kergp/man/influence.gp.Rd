\name{influence.gp}
\alias{influence.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
%% \encoding{UTF8} If ecole des Mines.. 
\title{
  Diagnostics for a Gaussian Process Model, Based on Leave-One-Out
}
\description{
  
  Cross Validation by leave-one-out for a \code{gp} object.

}
\usage{
\S3method{influence}{gp}(model, type = "UK", trend.reestim = TRUE, ...)
}
\arguments{
  \item{model}{

    An object of class \code{"gp"}.

  }
  \item{type}{

    Character string corresponding to the GP "kriging" family, to be
    chosen between simple kriging (\code{"SK"}), or universal kriging
    (\code{"UK"}).

  }
  \item{trend.reestim}{

    Should the trend be re-estimated when removing an observation?
    Default to \code{TRUE}.

  }
  \item{...}{

    Not used.

  }
}

\details{
  
  Leave-one-out (LOO) consists in computing the prediction at a design
  point when the corresponding observation is removed from the learning
  set (and this, for all design points). A quick version of LOO based on
  Dubrule's formula is also implemented; It is limited to 2 cases:
  \itemize{
    \item{\code{(type == "SK") & !trend.reestim} and}
    \item{\code{(type == "UK") & trend.reestim}.}
    }
}
  
\section{Warning}{
  
  Only trend parameters are re-estimated when removing one
  observation. When the number \eqn{n} of observations is small, the
  re-estimated values can be far away from those obtained with the
  entire learning set.
  
}

\value{

  A list composed of the following elements, where \emph{n} is the total
  number of observations.
  
  \item{mean }{

    Vector of length \emph{n}. The \eqn{i}-th element is the kriging
    mean (including the trend) at the \eqn{i}th observation number when
    removing it from the learning set.

  }
  \item{sd }{

    Vector of length \emph{n}. The \eqn{i}-th element is the kriging
    standard deviation at the \eqn{i}-th observation number when
    removing it from the learning set.
    
  }
}

\references{ 
  F. Bachoc (2013), "Cross Validation and Maximum Likelihood estimations of
  hyper-parameters of Gaussian processes with model
  misspecification". \emph{Computational Statistics and Data Analysis},
  \bold{66}, 55-69
  \href{http://www.sciencedirect.com/science/article/pii/S0167947313001187}{link}
  %% \url{http://www.lpma.math.upmc.fr/pageperso/bachoc/publications.html}
  
  N.A.C. Cressie (1993), \emph{Statistics for spatial data}. Wiley series
  in probability and mathematical statistics.
  
  O. Dubrule (1983), "Cross validation of Kriging in a unique
  neighborhood". \emph{Mathematical Geology}, \bold{15}, 687-699.

  J.D. Martin and T.W. Simpson (2005), "Use of kriging models to
  approximate deterministic computer models". \emph{AIAA Journal},
  \bold{43} no. 4, 853-863.
  
  M. Schonlau (1997), \emph{Computer experiments and global optimization}.
  Ph.D. thesis, University of Waterloo.
	
}

\author{

  O. Roustant, D. Ginsbourger.

}

\seealso{ \code{\link{predict.gp}},  \code{\link{plot.gp}} }

\keyword{models}
