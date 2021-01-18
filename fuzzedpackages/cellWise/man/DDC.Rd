\name{DDC}

\alias{DDC}

\title{
Detect Deviating Cells
}

\description{
This function aims to detect cellwise outliers in the data. These are entries in the data matrix which are substantially higher or lower than what could be expected based on the other cells in its column as well as the other cells in its row, taking the relations between the columns into account. Note that this function first calls \code{\link{checkDataSet}} and analyzes the remaining cleaned data.
}

\usage{
DDC(X, DDCpars = list())
}

\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or a data frame.
  }
  \item{DDCpars}{ A list of available options:
    \itemize{
      \item{ \code{fracNA} \cr
      Only consider columns and rows with fewer NAs (missing
      values) than this fraction (percentage). Defaults to \eqn{0.5}.}
      \item{ \code{numDiscrete} \cr
      A column that takes on \code{numDiscrete} or fewer values will
      be considered discrete and not used in the analysis. Defaults to \eqn{3}.}
      \item{ \code{precScale} \cr
      Only consider columns whose scale is larger than \code{precScale}.
      Here scale is measured by the median absolute deviation. Defaults to \eqn{1e-12}.}
      \item{ \code{cleanNAfirst} \cr
      If \code{"columns"}, first columns then rows are checked for NAs.
      If \code{"rows"}, first rows then columns are checked for NAs.
      \code{"automatic"} checks columns first if \eqn{d \geq 5n} and rows first otherwise.
      Defaults to \code{"automatic"}.
      }
      \item{ \code{tolProb} \cr
      Tolerance probability, with default \eqn{0.99}, which
      determines the cutoff values for flagging outliers in
      several steps of the algorithm.}
      \item{ \code{corrlim} \cr
      When trying to estimate \eqn{z_{ij}} from other variables \eqn{h}, we 
      will only use variables \eqn{h} with \eqn{|\rho_{j,h}| \ge corrlim}.
      Variables \eqn{j} without any correlated variables \eqn{h} satisfying 
      this are considered standalone, and treated on their own. Defaults to \eqn{0.5}.}
      \item{\code{combinRule} \cr
      The operation to combine estimates of \eqn{z_{ij}} coming from
      other variables \eqn{h}: can be \code{"mean"}, \code{"median"},
      \code{"wmean"} (weighted mean) or \code{"wmedian"} (weighted median).
      Defaults to \code{wmean}.}
      \item{ \code{returnBigXimp} \cr
      If TRUE, the imputed data matrix \code{Ximp} in the output
      will include the rows and columns that were not
      part of the analysis (and can still contain NAs). Defaults to \code{FALSE}.}
      \item{ \code{silent} \cr
      If \code{TRUE}, statements tracking the algorithm's progress will not be printed. Defaults to \code{FALSE}.}
       \item{ \code{nLocScale} \cr
    When estimating location or scale from more than \code{nLocScale} data values, the computation is based on a random sample of size \code{nLocScale} to save time. When  \code{nLocScale = 0} all values are used. Defaults to 25000.}
           \item{ \code{fastDDC} \cr
      Whether to use the fastDDC option or not. The fastDDC algorithm uses approximations
to allow to deal with high dimensions. Defaults to \code{TRUE} for \eqn{d > 750} and \code{FALSE} otherwise.}
      \item{ \code{standType} \cr
      The location and scale estimators used for robust standardization. Should be one of \code{"1stepM"}, \code{"mcd"} or \code{"wrap"}. See \code{\link{estLocScale}} for more info. Only used when \code{fastDDC = FALSE}. Defaults to \code{"1stepM"}.}
      \item{ \code{corrType} \cr
      The correlation estimator used to find the neighboring variables. Must be one of \code{"wrap"} (wrapping correlation), \code{"rank"} (Spearman correlation) or \code{"gkwls"} (Gnanadesikan-Kettenring correlation followed by weighting). Only used when \code{fastDDC} \code{ = FALSE}. Defaults to \code{"gkwls"}.}
            \item{ \code{transFun} \cr
      The transformation function used to compute the robust correlations when \code{fastDDC = TRUE}. Can be \code{"wrap"} or \code{"rank"}. Defaults to \code{"wrap"}.}
      \item{ \code{nbngbrs} \cr
     When \code{fastDDC = TRUE}, each column is predicted from at most \code{nbngbrs} columns correlated to it.
     Defaults to 100.}
      }
  }
}

\value{ A list with components: \cr
  \itemize{
    \item{\code{DDCpars} \cr
    The list of options used.
    }
    \item{\code{colInAnalysis} \cr
    The column indices of the columns used in the analysis.
    }
    \item{\code{rowInAnalysis} \cr
    The row indices of the rows used in the analysis.
    }
    \item{\code{namesNotNumeric} \cr
    The names of the variables which are not numeric.      
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
    Cleaned data after \code{checkDataSet}.           
    }
    \item{\code{locX} \cr
    Estimated location of \code{X}.           
    }
    \item{\code{scaleX} \cr
    Estimated scales of \code{X}.           
    }
    \item{\code{Z} \cr
    Standardized \code{remX}.           
    }
    \item{\code{nbngbrs} \cr
    Number of neighbors used in estimation.          
    }
    \item{\code{ngbrs} \cr
    Indicates neighbors of each column, i.e. the columns most correlated with it.   
    }
    \item{\code{robcors} \cr
    Robust correlations.           
    }
    \item{\code{robslopes} \cr
    Robust slopes.          
    }
    \item{\code{deshrinkage} \cr
    The deshrinkage factor used for every connected (i.e. non-standalone) column of \code{X}.           
    }              
    \item{\code{Xest} \cr
    Predicted \code{X}.           
    }       
     \item{\code{scalestres} \cr
    Scale estimate of the residuals \code{X - Xest}.
    }       
    \item{\code{stdResid} \cr
    Residuals of orginal \code{X} minus the estimated \code{Xest}, standardized by column.
    }
    \item{\code{indcells} \cr
    Indices of the cells which were flagged in the analysis.
    }
    \item{\code{Ti} \cr
    Outlyingness (test) value of each row.         
    }
    \item{\code{medTi} \cr
    Median of the Ti values.      
    }
    \item{\code{madTi} \cr
    Mad of the Ti values.       
    }
    \item{\code{indrows} \cr
    Indices of the rows which were flagged in the analysis.
    }
    \item{\code{indNAs} \cr
    Indices of all NA cells.
    } 
    \item{\code{indall} \cr
    Indices of all cells which were flagged in the analysis plus all cells in flagged rows plus the indices of the NA cells.
    }         
    \item{\code{Ximp} \cr
    Imputed \code{X}.           
    }
  }
}

\references{
Rousseeuw, P.J., Van den Bossche W. (2018). Detecting Deviating Data Cells. \emph{Technometrics}, \bold{60}(2), 135-145. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2017.1340909}{(link to open access pdf)}

Raymaekers, J., Rousseeuw P.J. (2019). Fast robust correlation for high dimensional data. \emph{Technometrics}, published online. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2019.1677270}{(link to open access pdf)}
}

\author{
 Raymaekers J., Rousseeuw P.J., Van den Bossche W. 
}


\seealso{
\code{\link{checkDataSet}},\code{\link{cellMap}}
}

\examples{
library(MASS); set.seed(12345)
n <- 50; d <- 20
A <- matrix(0.9, d, d); diag(A) = 1
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 50, FALSE)] <- NA
x[sample(1:(n * d), 50, FALSE)] <- 10
x[sample(1:(n * d), 50, FALSE)] <- -10
x <- cbind(1:n, x)
DDCx <- DDC(x)
cellMap(DDCx$remX, DDCx$stdResid,
columnlabels = 1:d, rowlabels = 1:n)

# For more examples, we refer to the vignette:
vignette("DDC_examples")
}
