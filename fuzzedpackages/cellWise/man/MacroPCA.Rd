\name{MacroPCA}

\alias{MacroPCA}

\title{
MacroPCA
}

\description{This function performs the MacroPCA algorithm, which can deal with Missing values and Cellwise
and Rowwise Outliers. Note that this function first calls  \code{\link{checkDataSet}} and analyzes the remaining cleaned data.}

\usage{
MacroPCA(X, k = 0, MacroPCApars = NULL)
}
\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or a data frame.
  }
  \item{k}{\code{k} is the desired number of principal components.
  If \code{k = 0} or \code{k = NULL}, the algorithm will compute the percentage
  of explained variability for \code{k} upto \code{kmax} and show a scree plot,
  and suggest to choose a value of k such that the cumulative percentage of
  explained variability is at least 80 \%.}
   
  \item{MacroPCApars}{A list of available options detailed below. If MacroPCApars = NULL the defaults below are used.
    \itemize{
      \item{\code{DDCpars} \cr A list with parameters for the first step of the MacroPCA
      algorithm (for the complete list see the function
      \code{\link{DDC}}). Default is \code{NULL}.}
   
      \item{\code{kmax} \cr The maximal number of principal components to compute. Default
       is \code{kmax = 10}. If \code{k} is provided kmax does not need to be specified,
       unless \code{k} is larger than 10 in which case you need to set \code{kmax}
       high enough.}
      
      \item{\code{alpha} \cr This is the coverage, i.e. the fraction of rows the algorithm
      should give full weight. Alpha should be between 0.50 and 1, the default is
      0.50.}
      
      \item{\code{scale} \cr A value indicating whether and how the original variables should
      be scaled. If \code{scale = FALSE} or \code{scale = NULL} no scaling is
      performed (and a vector of 1s is returned in the \code{$scaleX slot}).
      If \code{scale = TRUE} (default) the data are scaled by a 1-step M-estimator of scale with the Tukey biweight weight function to have a robust scale of 1.
      Alternatively scale can be a vector of length
      equal to the number of columns of \code{x}. The resulting scale estimates are
      returned in the \code{$scaleX} slot of the MacroPCA output.}
      
      \item{\code{maxdir} \cr The maximal number of random directions to use for computing the
      outlyingness of the data points. Default is \code{maxdir = 250}. If the number
      \eqn{n} of observations is small all \eqn{n * (n - 1) / 2} pairs of
      observations are used.}
      
      \item{\code{distprob} \cr The quantile determining the cutoff values
      for orthogonal and score distances. Default is 0.99.}
      
      \item{ \code{silent} \cr
      If TRUE, statements tracking the algorithm's progress will not be printed. Defaults to \code{FALSE}.}
      
      \item{\code{maxiter} \cr Maximum number of iterations. Default is 20.}
      
      \item{\code{tol} \cr Tolerance for iterations. Default is 0.005.}
      \item{\code{bigOutput} \cr whether to compute and return NAimp, Cellimp and Fullimp. Defaults to \code{TRUE}.}
    }
  }
}

\value{
  A list with components: \cr
  \item{MacroPCApars}{the options used in the call.}
  \item{remX}{Cleaned data after \code{checkDataSet}.}
  \item{DDC}{results of the first step of MacroPCA. These are needed to run
  MacroPCApredict on new data.}
  \item{scaleX}{the scales of the columns of \code{X}.}
  \item{k}{the number of principal components.}
  \item{loadings}{the columns are the \code{k} loading vectors.}
  \item{eigenvalues}{the \code{k} eigenvalues.}
  \item{center}{vector with the fitted center.}
  \item{alpha}{\code{alpha} from the input.}
  \item{h}{\code{h} (computed from \code{alpha}).}
  \item{It}{number of iteration steps.}
  \item{diff}{convergence criterion.}
  \item{X.NAimp}{data with all \code{NA}'s imputed by \code{MacroPCA}.}
  \item{scores}{scores of \code{X.NAimp}.}
  \item{OD}{orthogonal distances of the rows of \code{X.NAimp}.}
  \item{cutoffOD}{cutoff value for the OD.}
  \item{SD}{score distances of the rows of \code{X.NAimp}.}
  \item{cutoffSD}{cutoff value for the SD.}
  \item{indrows}{row numbers of rowwise outliers.}
  \item{residScale}{scale of the residuals.}
  \item{stdResid}{standardized residuals. Note that these are \code{NA}
  for all missing values of \code{X}.}
  \item{indcells}{indices of cellwise outliers.}
  \item{NAimp}{various results for the NA-imputed data.}
  \item{Cellimp}{various results for the cell-imputed data.}
  \item{Fullimp}{various result for the fully imputed data.}
}



\references{
Hubert, M., Rousseeuw, P.J., Van den Bossche W. (2019). MacroPCA: An all-in-one PCA method allowing for missing values as well as cellwise and rowwise outliers. \emph{Technometrics}, \bold{61}(4), 459-473. \href{https://www.tandfonline.com/doi/full/10.1080/00401706.2018.1562989}{(link to open access pdf)}
}

\author{
 Rousseeuw P.J., Van den Bossche W. 
}


\seealso{
\code{\link{checkDataSet}}, \code{\link{cellMap}},
\code{\link{DDC}}
}

\examples{
library(MASS) 
set.seed(12345) 
n <- 50; d <- 10
A <- matrix(0.9, d, d); diag(A) = 1
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 50, FALSE)] <- NA
x[sample(1:(n * d), 50, FALSE)] <- 10
x <- cbind(1:n, x)
MacroPCA.out <- MacroPCA(x, 2)
cellMap(MacroPCA.out$remX, MacroPCA.out$stdResid,
columnlabels = 1:d, rowlabels = 1:n)
# For more examples, we refer to the vignette:
vignette("MacroPCA_examples")
}