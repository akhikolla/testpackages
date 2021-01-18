\name{compute_A_B_G_D_and_simulate_mu_Lambda_q0_sameSigma}
\alias{compute_A_B_G_D_and_simulate_mu_Lambda_q0_sameSigma}
\title{
	Computation and simulations for \eqn{q = 0}.
}
\description{
	This function simulates \eqn{\mu}.
}
\usage{
	compute_A_B_G_D_and_simulate_mu_Lambda_q0_sameSigma(SigmaINV, 
		suff_statistics, K, priorConst1, T_INV, v_r)
}
\arguments{
  \item{SigmaINV}{
	Precision matrix \eqn{\Sigma^{-1}}
}
  \item{suff_statistics}{
	Sufficient statistics
}
  \item{K}{
	Number of overfitting mixture components
}
  \item{priorConst1}{
	Prior constant: \eqn{T^{-1}\xi}
}
  \item{T_INV}{
	Prior parameter: \eqn{T^{-1}\xi}
}
  \item{v_r}{
	This vector is used to set to zero the upper right \eqn{(q-1)\times(q-1)} diagonal block of factor loadings for identifiability purposes.
}
}
\value{
A list containing a draw from the conditional distributions of \eqn{\mu}:
	\item{mu}{\eqn{K\times p} array  (marginal mean per component)}

}
\author{
Panagiotis Papastamoulis
}


