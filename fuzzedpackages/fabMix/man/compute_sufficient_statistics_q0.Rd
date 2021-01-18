\name{compute_sufficient_statistics_q0}
\alias{compute_sufficient_statistics_q0}
\title{
Compute sufficient statistics for \eqn{q = 0}
}
\description{
Compute sufficient statistics given \eqn{z}.
}
\usage{
compute_sufficient_statistics_q0(z, K, x_data)
}
\arguments{
  \item{z}{
	Allocation vector
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	Data
}
}
\value{
A list with six entries of sufficient statistics.
	\item{cluster_size}{ Integer vector of length \eqn{K} }
	\item{sx}{ \eqn{K\times p} array }
	\item{sy}{ Not used here }
	\item{sxx}{ Not used }
	\item{syy}{ Not used here }
	\item{sxy}{ Not used here }

}
\author{
Panagiotis Papastamoulis
}
\examples{
        data(waveDataset1500)
        x_data <- as.matrix(waveDataset1500[ 1:20, -1]) # data
        z <-  waveDataset1500[ 1:20, 1] # class
        p <- dim(x_data)[2]
        n <- dim(x_data)[1]
        q <- 2
        K <- length(table(z))           # 3 classes
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_q0(
	 z = z, K = K, x_data = x_data)
}

