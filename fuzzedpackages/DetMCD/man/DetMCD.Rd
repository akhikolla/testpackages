\name{DetMCD}
\alias{DetMCD}
\title{Robust and Deterministic Location and Scatter Estimation via DetMCD}
\description{
     Computes a robust and deterministic multivariate location and scatter 
	estimate with a high breakdown point, using the \code{DetMCD} (Deterministic 
    Minimum Covariance Determinant) algorithm.
}
\usage{
  DetMCD(X,h=NULL,alpha=0.75,scale_est="Auto",tol=1e-07)
}
\arguments{
  \item{X}{a numeric matrix or data frame.
 Missing values (NaN's) and infinite values (Inf's) are allowed: observations (rows)
 with missing or infinite values will automatically be excluded from the computations.}
  \item{alpha}{Ignored if \code{h!=NULL}.  (Possibly vector of) numeric parameter controlling the size of the subsets over
          which the determinant is minimized, i.e., \code{alpha*n}
          observations are used for computing the determinant.  Allowed
          values are between 0.5 and 1 and the default is 0.75.}
  \item{h}{numeric integer parameter controlling the size of the subsets over
          which the determinant is minimized, i.e., \code{h}
          observations are used for computing the determinant.  Allowed
          values are between \code{[(n+p+1)/2]} and \code{n} and the default is NULL.}
  \item{scale_est}{a character string specifying the
  variance functional.  Possible values are "qn", "tau" and 'Auto".
Default value \code{"Auto"} is to use the Qn
              estimator for data with less than 1000 observations, and to use the
              tau-scale for data sets with more observations. But one
              can also always use the Qn estimator \code{"qn"}
              or the tau scale \code{"tau"}.}

  \item{tol}{a small positive numeric value to be
  used for determining numerical 0.}
}
\value{
  A list with components:


     \item{raw.center}{ The raw MCD location of the data.}
        \item{raw.cov}{The raw MCD covariance matrix (multiplied by a
                  consistency factor).}
  \item{crit}{The determinant of the raw MCD covariance matrix.}
        \item{raw.rd}{The robust distance of each observation to the raw MCD center, relative to 
                  the raw MCD scatter estimate. }
         \item{raw.wt}{ Weights based on the estimated raw covariance matrix 'raw.cov' and
                  the estimated raw location 'raw.center' of the data. These weights determine
                  which observations are used to compute the final MCD estimates.}

       \item{center}{The robust location of the data, obtained after
                    reweighting.}
          \item{cov}{The robust covariance matrix, obtained after
                    reweighting.}
          \item{h}{The number of observations that have determined the MCD estimator,
                    i.e. the value of h.}
                    \item{which.one}{The identifier 
                        of the initial shape estimate which led to the
                        optimal result.}  
                    \item{best}{ The subset of h points whose covariance matrix has minimal determinant.}
\item{weights}{ The finale vector of weights.}
           \item{rd}{ The robust distance of each observation to the final,
                   reweighted MCD center of the data, relative to the
                    reweighted MCD scatter of the data.  These distances allow
                    us to easily identify the outliers.}
          \item{rew.md}{ The Mahalanobis distance of each observation (distance from the classical
                    center of the data, relative to the classical shape
                    of the data).}
            \item{X}{ Same as the X in the call to DetMCD,
                    without rows containing missing or infinite values.}
	\item{alpha}{The vector of values of alpha used in the algorithm.}
	\item{scale_est}{The vector of scale estimators used in the estimates (one of \code{tau2} or \code{qn}.}
}
\details{
   DetMCD computes the MCD estimator of a multivariate data set in a deterministic way.
 This estimator is given by the subset of h observations with smallest
 covariance determinant.  The MCD location estimate is then the mean of those h points,
 and the MCD scatter estimate is their covariance matrix.  The default value
 of h is roughly 0.75n (where n is the total number of observations), but the
 user may choose each value between n/2 and n. Based on the raw estimates,
 weights are assigned to the observations such that outliers get zero weight.
 The reweighted MCD estimator is then given by the mean and covariance matrix
 of the cases with non-zero weight. 

 To compute the MCD estimator, six initial robust h-subsets are
 constructed based on robust transformations of variables or robust and
 fast-to-compute estimators of multivariate location and shape. Then
 C-steps are applied on these h-subsets until convergence. Note that the
 resulting algorithm is not fully affine equivariant, but it is often
 faster than the FAST-MCD algorithm which is affine equivariant. 
 Note that this function can not handle exact fit
 situations: if the raw covariance matrix is singular, the program is
 stopped. In that case, it is recommended to apply the FastMCD function.

 
 The MCD method is intended for continuous variables, and assumes that
 the number of observations n is at least 5 times the number of variables p.
 If p is too large relative to n, it would be better to first reduce
 p by variable selection or robust principal components (see the functions
 PcaHubert).
}
\references{
   Hubert, M., Rousseeuw, P.J. and Verdonck, T. (2012),
   "A deterministic algorithm for robust location and scatter", Journal of
   Computational and Graphical Statistics, Volume 21, Number 3, Pages 618--637.

   Verboven, S., Hubert, M. (2010). Matlab library LIBRA, Wiley Interdisciplinary Reviews: Computational Statistics, 2, 509--515.
}
\examples{
## generate data
set.seed(1234)  # for reproducibility
alpha<-0.5
n<-101
p<-5
#generate correlated data
D<-diag(rchisq(p,df=1))
W<-matrix(0.9,p,p)
diag(W)<-1
W<-D%*%W%*%t(D)
V<-chol(W)
x<-matrix(rnorm(n*p),nc=p)
x<-scale(x)%*%V


result<-DetMCD(x,scale_est="tau",alpha=alpha)
plot(result, which = "dd")

#compare to robustbase:
result<-DetMCD(x,scale_est="qn",alpha=alpha)
resultsRR<-covMcd(x,nsamp='deterministic',scalefn=qn,alpha=alpha)
#should be the same:
result$crit
resultsRR$crit


#Example with several values of alpha:
alphas<-seq(0.5,1,l=6)
results<-DetMCD(x,scale_est="qn",alpha=alphas)
plot(results, h.val = 2, which = "dd")
}
\author{
  Vakili Kaveh (includes section of the help file from the LIBRA implementation).
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
