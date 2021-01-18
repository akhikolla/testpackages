\name{checkDataSet}
\alias{checkDataSet}

\title{
Clean the dataset
}

\description{
This function checks the dataset X, and sets aside certain
columns and rows that do not satisfy the conditions.
It is used by the \code{\link{DDC}} and \code{\link{MacroPCA}} functions but can be used by itself, to clean a dataset for a different type of analysis.
}
\usage{
checkDataSet(X, fracNA = 0.5, numDiscrete = 3, precScale = 1e-12, silent = FALSE,
cleanNAfirst = "automatic")
}

\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or data frame.
}
  \item{fracNA}{
   Only retain columns and rows with fewer NAs than this fraction.
   Defaults to \eqn{0.5}.
}
  \item{numDiscrete}{
  A column that takes on numDiscrete or fewer values
  will be considered discrete and not retained in the cleaned data.
  Defaults to \eqn{3}.
}
  \item{precScale}{
  Only consider columns whose scale is larger than precScale.
  Here scale is measured by the median absolute deviation.
  Defaults to \eqn{1e-12}.
}
 \item{silent}{
  Whether or not the function progress messages should be printed.
  Defaults to \code{FALSE}.
}
 \item{cleanNAfirst}{
  If \code{"columns"}, first columns then rows are checked for NAs.
  If \code{"rows"}, first rows then columns are checked for NAs.
  \code{"automatic"} checks columns first if \eqn{d \geq 5n} and rows first otherwise.
  Defaults to \code{"automatic"}.
}
}

\value{
A list with components: \cr
  \itemize{
    \item{\code{colInAnalysis} \cr
    Column indices of the columns used in the analysis.
    }
    \item{\code{rowInAnalysis}} \cr
    Row indices of the rows used in the analysis.
    \item{\code{namesNotNumeric} \cr
    Names of the variables which are not numeric.       
    }
    \item{\code{namesCaseNumber} \cr
    The name of the variable(s) which contained the case numbers and was therefore removed.
    }
    \item{\code{namesNAcol} \cr
    Names of the columns left out due to too many \code{NA}'s.
    }
    \item{\code{namesNArow} \cr
    Names of the rows left out due to too many \code{NA}'s.
    }
    \item{\code{namesDiscrete} \cr
    Names of the discrete variables.           
    }
    \item{\code{namesZeroScale} \cr
    Names of the variables with zero scale.           
    }
    \item{\code{remX} \cr
    Remaining (cleaned) data after checkDataSet.           
    }
  }
}
\references{
Rousseeuw, P.J., Van den Bossche W. (2018). Detecting Deviating Data Cells. \emph{Technometrics}, \bold{60}(2), 135-145. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2017.1340909}{(link to open access pdf)}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}

\seealso{
\code{\link{DDC}}, \code{\link{MacroPCA}}, \code{\link{transfo}}, \code{\link{wrap}}
}

\examples{
library(MASS) 
set.seed(12345) 
n <- 100; d = 10
A <- matrix(0.9, d, d); diag(A) = 1
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 100, FALSE)] <- NA
x <- cbind(1:n, x)
checkedx <- checkDataSet(x)

# For more examples, we refer to the vignette:
vignette("DDC_examples")
}
