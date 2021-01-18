\name{DDCpredict}

\alias{DDCpredict}

\title{
DDCpredict
}

\description{Based on a \code{\link{DDC}} fit on an initial (training) data set \code{X}, this function
analyzes a new (test) data set \code{Xnew}.}

\usage{
DDCpredict(Xnew, InitialDDC, DDCpars = NULL)
}

\arguments{

\item{Xnew}{The new data (test data), which must be a matrix or a data frame.
It must always be provided.}

\item{InitialDDC}{The output of the \code{\link{DDC}} function on the initial (training)
dataset. Must be provided.}

\item{DDCpars}{The input options to be used for the prediction.
By default the options of InitialDDC are used.}
  
}

\value{
  A list with components: \cr
  \item{DDCpars}{the options used in the call, see \code{\link{DDC}}.}
  \item{locX}{the locations of the columns, from \code{InitialDDC}.}
  \item{scaleX}{the scales of the columns, from \code{InitialDDC}.}
  \item{Z}{\code{Xnew} standardized by \code{locX} and \code{scaleX}.}
  \item{nbngbrs}{predictions use a combination of \code{nbngbrs} columns.}
  \item{ngbrs}{for each column, the list of its neighbors, from \code{InitialDDC}.}
  \item{robcors}{for each column, the correlations with its neighbors, from \code{InitialDDC}.}
  \item{robslopes}{slopes to predict each column by its neighbors, from \code{InitialDDC}.}
  \item{deshrinkage}{for each connected column, its deshrinkage factor used in \code{InitialDDC}.}
  \item{Xest}{predicted values for every cell of \code{Xnew}.}
  \item{scalestres}{scale estimate of the residuals (\code{Xnew} - \code{Xest}), from \code{InitialDDC}.}
  \item{stdResid}{columnwise standardized residuals of \code{Xnew}.}
  \item{indcells}{positions of cellwise outliers in \code{Xnew}.}
  \item{Ti}{outlyingness of rows in \code{Xnew}.}
  \item{medTi}{median of the \code{Ti} in \code{InitialDDC}.}
  \item{madTi}{mad of the \code{Ti} in \code{InitialDDC}.}
  \item{indrows}{row numbers of the outlying rows in \code{Xnew}.}
  \item{indNAs}{positions of the \code{NA}'s in \code{Xnew}.}
  \item{indall}{positions of \code{NA}'s and outlying cells in \code{Xnew}.}
  \item{Ximp}{\code{Xnew} where all cells in indall are imputed by their prediction.}
  
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
n <- 100; d <- 10
A <- matrix(0.9, d, d); diag(A) = 1
x <- mvrnorm(n, rep(0,d), A)
x[sample(1:(n * d), 50, FALSE)] <- NA
x[sample(1:(n * d), 50, FALSE)] <- 10
x <- cbind(1:n, x)
DDCx <- DDC(x)
xnew <- mvrnorm(50, rep(0,d), A)
xnew[sample(1:(50 * d), 50, FALSE)] <- 10
predict.out <- DDCpredict(xnew, DDCx)
cellMap(xnew, predict.out$stdResid,
columnlabels = 1:d, rowlabels = 1:50)

# For more examples, we refer to the vignette:
vignette("DDC_examples")
}
