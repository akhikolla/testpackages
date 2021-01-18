\name{plot.clust}
\alias{plot.clust}
\title{Plots clustering result}
\usage{
  \method{plot}{clust} (x,
    type = c("histo", "scat", "sil"), adjust = TRUE, ...)
}
\arguments{
  \item{x}{The clustering results obtained from
  \code{\link{ClusProc}}.}

  \item{type}{Factor. For specifying the plot type. It must
  be one of 'histo', 'scat' and 'sil'. If it is 'histo',
  the histogram is obtained with the first PC score of the
  intensity measurement. For 'scat', the first PC score of
  the intensity measurement is plotted against the mean of
  the intensity measurement. For 'sil', the silhouette
  score is plotted. See details.}

  \item{adjust}{Logicals. If TRUE (default), the
  silhouette-adjusted clustering result will be used. If
  FALSE, the initial clustering result will be used. See
  details in \code{\link{ClusProc}}.}

  \item{...}{Usual arguments passed to the qplot function.}
}
\description{
  Makes formatted plots from the clustering result returned
  from \code{\link{ClusProc}}.
}
\details{
  \itemize{ \item{type}{We provide three types of plots:
  'hist', 'scat' and 'sil'. The first two plots are used to
  visually check the performance of clustering. Different
  clusters are represented by using different colors. The
  'sil' plot is the the overview of the silhouette value
  for all the individuals, the silhouettes of the different
  clusters are printed below each other. The higher
  silhouettes value means the better performance.} }
}
\examples{
# Fit the data under the given clustering numbers
clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
plot(clus.fit,type='histo')
}
\author{
  Meiling Liu
}

