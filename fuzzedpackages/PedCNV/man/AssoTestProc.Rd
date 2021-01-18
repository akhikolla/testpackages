\name{AssoTestProc}
\alias{AssoTestProc}
\title{CNV association test procedure}
\usage{
  AssoTestProc(signal, fam, envirX, phi, N,
    varSelection = c("PC1", "RAW", "PC.9", "MEAN"),
    H0 = TRUE, threshold = 1e-05, itermax = 8,
    thresEM = 0.005, thresAI = 1e-05)
}
\arguments{
  \item{signal}{The matrix of intensity measurements. The
  row names must be consistent with the Individual ID in
  fam file.}

  \item{fam}{The FAM file which follows the format defined
  in PLINK.}

  \item{envirX}{The matrix of environmental variables. The
  intercept is automatically included and it does not need
  to be in this matrix.}

  \item{phi}{The correlation matrix between individuals. It
  can be built with the kinship coefficient or the
  estimated correlation matrix with SNP data. Free software
  that builds this matrix is available, and one of them can
  be downloaded at \url{http://biostat.ac.kr/fqls/} The
  default is an identity matrix and it is for independent
  samples.}

  \item{N}{Number of clusters one wants to fit to the data.
  N needs to be larger than 1 and if it is 1, error will be
  returned. It can be estimated with the function
  \code{\link{ClusProc}}.}

  \item{varSelection}{Factor. For specifying how to handle
  the intensity values. It must take value on 'RAW',
  'PC.9', 'PC1'and 'MEAN'. If the value is 'RAW', then the
  raw intensity value will be used. If it is 'PC.9', then
  the first several PCA scores which account for 90\% of
  all the variance will be used. If the value is 'PC1',
  then the first PCA scores will be used. If the value is
  'MEAN', the mean of all the probes will be used. The
  default method is 'PC1'.}

  \item{H0}{Logicals. If it is TRUE (the default), all
  parameters are estimated under the assumption that there
  is no genetic association between CNV and phenotypes. If
  it is FALSE, parameters are estimated under the null or
  alternative hypothesis.}

  \item{threshold}{Optional number of convergence
  threshold. The iteration stops if the absolute difference
  of log likelihood between successive iterations is less
  than it. The default threshold 1e-05 will be used if it's
  missing.}

  \item{itermax}{Optional. The iteration stops if the times
  of iteration is large than this value. The default number
  8 will be used if it's missing.}

  \item{thresEM}{Optional number of convergence threshold
  in the EM (expectation-maximization method) procedure.
  The default threshold 0.005 will be used if it's
  missing.}

  \item{thresAI}{Optional number of convergence threshold
  in the AI (average information method) procedure. The
  default threshold 1e-05 will be used if it's missing.}
}
\value{
  It returns object of class 'asso'. The result is obtained
  under the null hypothesis if H0 is TRUE, otherwise the
  result is obtained under null or alternative hypothesis.
  \item{para}{The parameter estimations for the best fit.}
  \item{clusRes}{The clustering assignment for each
  individual.}
}
\description{
  This function tests the association of CNV with
  continuous trait of interest. Two statistics are provided
  for different strategies with the intensity measurement.
}
\examples{
# Fit the data under the assumption that there are 3 clusters
fit.pc <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
}
\author{
  Meiling Liu, Sungho Won and Weicheng Zhu
}

