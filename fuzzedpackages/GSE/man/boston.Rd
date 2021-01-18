\name{boston}
\alias{boston}
\title{Boston Housing Data}
\usage{ data(boston) }
\description{Housing data for 506 census tracts of Boston from the 1970
  census. The dataframe \code{boston} contains the corrected data by Harrison and
  Rubinfeld (1979). The data was for a few minor errors and augmented with the 
  latitude and longitude of the observations.
  The original data can be found in the references below.
}
\format{The original data are 506 observations on 14 variables,
    \code{medv} being the target variable:
    \tabular{ll}{
      \code{cmedv} \tab corrected median value of owner-occupied homes in USD 1000's \cr
      \code{crim} \tab per capita crime rate by town \cr
      \code{indus} \tab proportion of non-retail business acres per town \cr
      \code{nox} \tab nitric oxides concentration (parts per 10 million) \cr
      \code{rm} \tab average number of rooms per dwelling \cr
      \code{age} \tab proportion of owner-occupied units built prior to 1940 \cr
      \code{dis} \tab weighted distances to five Boston employment centres \cr
      \code{rad} \tab index of accessibility to radial highways \cr
      \code{tax} \tab full-value property-tax rate per USD 10,000 \cr
      \code{ptratio} \tab pupil-teacher ratio by town \cr
      \code{b} \tab \eqn{1000(B - 0.63)^2} where \eqn{B} is the proportion of blacks by town\cr
      \code{lstat} \tab percentage of lower status of the population \cr
    }
}
\references{  
  Harrison, D. and Rubinfeld, D.L. (1978).
  Hedonic prices and the demand for clean air.
  \emph{Journal of Environmental Economics and Management}, \bold{5},
  81--102.

  Gilley, O.W., and R. Kelley Pace (1996). On the Harrison and Rubinfeld  
  Data. \emph{Journal of Environmental Economics and Management}, \bold{31},  
  403--405. [Provided corrections and examined censoring.]

  Newman, D.J. & Hettich, S. & Blake, C.L. & Merz, C.J. (1998).
  UCI Repository of machine learning databases
  [http://www.ics.uci.edu/~mlearn/MLRepository.html]. Irvine, CA:
  University of California, Department of Information and Computer
  Science.

  Pace, R. Kelley, and O.W. Gilley (1997). Using the Spatial Configuration of  
  the Data to Improve Estimation. \emph{Journal of the Real Estate Finance  
  and Economics}, \bold{14}, 333--340. [Added georeferencing and spatial
  estimation.]
}  
\source{
  The original data have been taken from the UCI Repository Of Machine Learning
  Databases at
  \itemize{
    \item \url{http://www.ics.uci.edu/~mlearn/MLRepository.html},
  }
  the corrected data have been taken from Statlib at
  \itemize{
    \item \url{http://lib.stat.cmu.edu/DASL/}
  }
  See Statlib and references there for details on the corrections.
  Both were converted to R format by Friedrich Leisch.
}
\keyword{datasets}

