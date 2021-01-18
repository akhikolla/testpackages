\name{shape.fd.outliers}
\alias{shape.fd.outliers}
\title{Functional Depth-Based Shape Outlier Detection}
\usage{
shape.fd.outliers(dataf, range = NULL, d = 101, q = 0.05,
  method = c("halfspace", "simplicial"), approx = 100, print = FALSE,
  plotpairs = FALSE, max.order = 3, exclude.out = TRUE,
  output = c("matrix", "list"), identifiers = NULL)
}
\arguments{
\item{dataf}{Functional dataset, represented by a \code{dataf} object of their arguments
and functional values. \code{n} stands for the number of functions.}

\item{range}{The common range of the domain where the fucntions \code{dataf} are observed.
Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
\code{dataf}.}

\item{d}{Grid size to which all the functional data are transformed. For depth computation, 
all functional observations are first transformed into vectors of their functional values of length \code{d}
corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
points are reconstructed using linear interpolation, and extrapolation.}

\item{q}{The quantile presenting a threshold for the first order outlier detection. Functions with first order integrated depth
smaller than the \code{q} quantile of this sample of depths are flagged as potential outliers. If set to \code{NULL}, the
the outliers are detected from the first order integrated depth after the log-transformation, as for higher order outliers.}

\item{method}{The depth that is used in the diagnostic plot. possible values are \code{halfspace} for 
the halfspace depth, or \code{simplicial} for the simplicial depth.}

\item{approx}{For the computation of the third order integrated depth,
the number of approximations used in the computation of the order extended depth. By default
this is set to \code{100}, meaning that \code{100}
trivariate points are randomly sampled in unit cube, and at these points the trivariate depths of the
corresponding functional values. May be set to \code{0} to compute the depth at all possible \code{d^3}
combinations of the points in the domain. This choice may result in very slow computation, see also \code{\link{depthf.fd1}}.}

\item{print}{If the rows of \code{X} are named, \code{print=TRUE} enables a graphical output when the names of the outlying curves
are displayed.}

\item{plotpairs}{If set to \code{TRUE}, the scatter plot of the computed depths for orders \code{1}, \code{2} and \code{3} is
is displayed. Here, the depths corresponding to the flagged outliers are plotted in colour.}

\item{max.order}{Maximal order of shape outlyingness to be computed, can be set to \code{1}, \code{2}, or \code{3}.}

\item{exclude.out}{Logical variable; exclude the detected lower order outliers in the flagging process? By default \code{TRUE}.}

\item{output}{Output method, can be set to \code{matrix} for a matrix with logical entries (\code{TRUE} for outliers), or \code{list} for 
a list of outliers.}

\item{identifiers}{A vector of names for the data observations. Facilitates identification of outlying functions.}
}
\value{
A matrix of logical values of size \code{n*4}, where \code{n} is the sample size. In the first three rows indicators of outlyingness
of the corresponding functions for orders \code{1}, \code{2} and \code{3} are given, in the fourth row the indicator of outlyingness
with respect to the comparison of the first, and third order depths is given. That is, the fist row corresponds to the first order outliers, 
the second row to the second order outliers, and the last two rows formally to the third order outliers. Please consult Nagy et al. (2016)
to interpret the notion of shape outlyingness.
}
\description{
Detects functional outliers of first three orders, based on the order extended integrated depth for functional data.
}
\details{
Using the procedure described in Nagy et al. (2016), the function uses the order extended integrated depths for functions, 
see \code{\link{depthf.fd1}} and \code{\link{shape.fd.analysis}}, to perform informal functional shape outlier detection. 
Outliers of the first order (horizontal shift outliers) are found as the functions with \code{q} \% of smallest (first order)
integrated depth values. Second and third order outliers (shape outliers) are found using the extension of the boxplot method
for depths as described in the paper Nagy et al. (2016).
}
\examples{
n = 30
dataf = dataf.population()$dataf[1:n]
shape.fd.outliers(dataf,print=TRUE,plotpairs=TRUE,
identifiers=unlist(dataf.population()$identifier)[1:n])

}
\references{
Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
Depth-based recognition of shape outlying functions. 
\emph{Journal of Computational and Graphical Statistics}, \bold{26} (4), 883--893.
}
\seealso{
\code{\link{depthf.fd1}}, \code{\link{shape.fd.analysis}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
\keyword{outlier}
