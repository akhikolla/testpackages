\name{lsmacof}
\alias{lsmacof}

\title{
3D reconstruction from Hi-C distances
}
\description{
A function for reconstructing the 3D configuration from pairwise calibrated Hi-C distances through local stress minimization (see details).  
}
\usage{
lsmacof(diss, infD, itermax = 10000, eps = 1e-06, init = NULL, k = 3,
verbose = FALSE, infW = NULL)
}
\arguments{
  \item{diss}{
A M*M Hi-C matrix,providing pairwise calibrated distances between M genomic bins. Zero off-diagonal values represent discarded/missing information. This matrix could be prepared with the function \code{\link{calibrate}}.
}
  \item{infD}{
A numeric value for missing distances (see details).
}
  \item{itermax}{
An integer value giving the maximal number of iterations for the SMACOF procedure (see details).Set to be 10000 by default.
}
  \item{eps}{
A numeric value giving the convergence parameter for the SMACOF procedure (see details). Set to be 1e-06 be default. 
}
  \item{init}{
A M*k numeric matrix, giving the initial configuration, for the SMACOF procedure, and set to NULL by default. If NULL, the classical MDS solution is used.
}
  \item{k}{
The number of dimensions for the output configuration. Set to 3 by default.
}
  \item{verbose}{
A Boolean indicating whether to print the stress at each iteration.Set to FALSE by default.
}
  \item{infW}{
A numeric value giving the weight for missing distances. Set to NULL by default. If NULL, this value is set to be 1/infD for infD > 1 and 0.05 otherwise (see details).
}
}
\details{
A calibrated distance Hi-C matrix H can be used as input for a 3D reconstruction algorithm. It is important to note that if H would give the true 3D Euclidean distances than the multi-dimensional scaling (MDS) solution (Torgerson, 1952) would recover the underlying 3D configuration, up to a rotation. However, due to Hi-C limitation we rely mostly on local calibrated distances. Denote delta(i,j) the calibrated distance between loci i and j, and d(i,j), their Euclidean distance, in the true underlying 3D configuration. Here, a zero delta[i,j] , for different loci i and j, represents a missing information that was discarded in the calibration step. Our goal is then to minimize the following function, usually termed stress in an MDS setting (Kruskal, 1964): SUM(i<j)[w(i,j)(delta(i,j)-d(i,j)(Y))^2] , where w(i,j) are the weights we assign according to the reliability of delta(i,j). Since we mostly rely on local information, we can use here a local stress function (Chen and Buja, 2009), where missing delta(i,j) are replaced with a constant dinf  (dinf >> known delta(i,j)) and w(i,j) take the value of 1/dinf for missing distances and 1 otherwise (for dinf equal or smaller than 1, weights of missing distances should be set to a small constant << 1). Since w(i,j) define an irreducible matrix, the stress minimization could be performed through Scaling by Majorizing a Complicated Function (SMACOF) (De Leeuw, 1977; De Leeuw, 1988; De Leeuw and Heiser 1977), a well-established strategy for this task, which guarantees convergence.
}
\value{
A M*k configuration matrix (k=3 by default) reconstructed from the given M*M Hi-C pairwise distances.

}
\references{
Main reference:\cr
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).\cr

References for Details section:\cr
Chen, L. and Buja, A. (2009) Local Multidimensional Scaling for Nonlinear Dimension Reduction, Graph Drawing, and Proximity Analysis. Journal of the American Statistical Association 104, 209{-}219.\cr
De Leeuw, J. (1977). Applications of convex analysis to multidimensional scaling. In Barra, J.R et al (Eds.), Recent developments in statistics, 133{-}145. Amsterdam, The Netherlands: North-Holland.\cr
De Leeuw, J. (1988). Convergence of the majorization method for multidimensional scaling,. Journal of Classification, 5, 163{-}180.\cr
De Leeuw, J. and Heiser, W.J. (1977). Convergence of correction-matrix algorithms for multidimensional scaling. In Lingoes, J.C., Roskam, E.E., and Borg, I. (Eds.),Geometric representations of relational data, 735{-}752. Ann Arbor, MI:Mathesis.\cr
Kruskal, J. B. (1964). Nonmetric multidimensional scaling: A numerical method. Psychometrika, 29, 115{-}129.\cr
Torgerson, W.S. (1952). Multidimensional scaling: I. Theory and method. Psychometrika, 17, 401{-}19.\cr
}
\author{
Yoli Shavit
}
\seealso{
\code{\link{prepareCalib}}\cr
\code{\link{calibHiC}}\cr
\code{\link{calibrate}}
}

\examples{
 data(calibHiC)
 data(match)
 data(conf)
 predConf = lsmacof(calibHiC, max(match$distances))
 
 # superimpose
 partialPS<-function(m1, m2)
 {
    # translate to origin
    tm1<-scale(m1, scale = FALSE)
    tm2<-scale(m2, scale = FALSE)
    A<-svd(t(tm2)\%*\%tm1)
    v<-A$u
    w<-A$v
    # update v a det(R) is positive
    k = ncol(m1)
    d = sign(det(t(w)\%*\%t(v)))
    v[,k] = v[,k]*-1*d
    R<- w\%*\%t(v)
    return(list(m1=tm1\%*\%R,m2=tm2))
  }

  res = partialPS(predConf, conf) 
  if (require(rgl))
  {
    plot3d(res$m2, type = "l", axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "")
    lines3d(res$m1, col = "red")
  }

}
