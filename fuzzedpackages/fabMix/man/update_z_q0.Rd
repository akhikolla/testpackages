\name{update_z_q0}
\alias{update_z_q0}
\title{
Gibbs sampling for \eqn{z} for \eqn{q=0}
}
\description{
Gibbs sampling for \eqn{z}
}
\usage{
update_z_q0(w, mu, SigmaINV, K, x_data)
}
\arguments{
  \item{w}{
	Mixture weights
}
  \item{mu}{
	Marginal means
}
  \item{SigmaINV}{
	Precision matrix per component
}
  \item{K}{
	Number of components
}
  \item{x_data}{
	Data
}
}
\value{
	Allocation vector
}
\author{
Panagiotis Papastamoulis
}

