\name{ClusProc}
\alias{ClusProc}
\title{CNV clustering Procedure}
\usage{
  ClusProc(signal, N = 2:6,
    varSelection = c("PC1", "RAW", "PC.9", "MEAN"),
    threshold = 1e-05, itermax = 8, adjust = TRUE,
    thresMAF = 0.01, scale = FALSE, thresSil = 0.01)
}
\arguments{
  \item{signal}{The matrix of intensity measurements. The
  row names must be consistent with the Individual ID in
  fam file.}

  \item{N}{Number of clusters one wants to fit to the data.
  N needs to be larger than 1 and if it is 1, error will be
  returned. The default value 2,3,...,6 will be used if it
  is missing.}

  \item{varSelection}{Factor. For specifying how to handle
  the intensity values. It must take value on 'RAW',
  'PC.9', 'PC1'and 'MEAN'. If the value is 'RAW', then the
  raw intensity value will be used. If it is 'PC.9', then
  the first several PCA scores which account for 90\% of
  all the variance will be used. If the value is 'PC1',
  then the first PCA scores will be used. If the value is
  'MEAN', the mean of all the probes will be used. The
  default method is 'PC1'.}

  \item{threshold}{Optional number of convergence
  threshold. The iteration stops if the absolute difference
  of log likelihood between successive iterations is less
  than it. The default threshold 1e-05 will be used if it's
  missing.}

  \item{itermax}{Optional. The iteration stops if the time
  of iteration is large than this value. The default number
  8 will be used if it's missing.}

  \item{adjust}{Logicals, If TRUE (default), the result
  will be adjusted by the silhouette score. See details.}

  \item{thresMAF}{The minor allele frequency threshold.}

  \item{thresSil}{The abandon threshold. The individual
  whose silhouette score is smaller than this value will be
  abandoned.}

  \item{scale}{Logicals. If TRUE, the signal will be scale
  by using sample mean and sample variance by columns
  before further data-processing.}
}
\value{
  It returns object of class 'clust'. 'clust' is a list
  containing following components: \item{clusNum}{The
  optimal number of clusters among give parameter {N}.}
  \item{silWidth}{Silhouette related results.}
}
\description{
  This function chooses the optimal number of clusters and
  provides the assignments of each individuals under the
  optimum clustering number.
}
\details{
  \itemize{ \item{adjust}{If adjust is TRUE, the result
  will be adjusted by the silhouette score in the following
  criterion. For each individual, the silhouette scores are
  calculated for each group. The individual will assigned
  forcefully to the group which maximize the silhouette
  scores. } }
}
\examples{
# Fit the data under the given clustering numbers
clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
}
\author{
  Meiling Liu
}

