\name{compute_A_B_G_D_and_simulate_mu_Lambda_CUU}
\alias{compute_A_B_G_D_and_simulate_mu_Lambda_CUU}
\title{
	Computation and simulations for CUU
}
\description{
	This function simulates \eqn{\mu} and \eqn{\Lambda} for the CUU model.
}
\usage{
	compute_A_B_G_D_and_simulate_mu_Lambda_CUU(SigmaINV, 
		suff_statistics, OmegaINV, K, priorConst1, T_INV, v_r)
}
\arguments{
  \item{SigmaINV}{
	Precision matrix \eqn{\Sigma^{-1}}
}
  \item{suff_statistics}{
	Sufficient statistics
}
  \item{OmegaINV}{
	Prior parameter: \eqn{\Omega^{-1}}
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
A list containing a draw from the conditional distributions of \eqn{\mu} and \eqn{\Lambda}:
	\item{Lambdas}{\eqn{K\times p\times q} array (factor loadings per component)}
	\item{mu}{\eqn{K\times p} array  (marginal mean per component)}

}
\author{
Panagiotis Papastamoulis
}
\examples{
	library('fabMix')
	data(waveDataset1500)
	x_data <- scale(as.matrix(waveDataset1500[ 1:20, -1])) # data
	z <-  waveDataset1500[ 1:20, 1] # class
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	q <- 2
	K <- length(table(z))           # 3 classes
	# give some arbitrary values to the parameters:
	set.seed(1)
	mu <- array( runif(K * p), dim = c(K,p) )
	y <- array(rnorm(n = q*n), dim = c(n,q))
	SigmaINV <- array(data = 0, dim = c(K,p,p) )
	for(k in 1:K){
		diag(SigmaINV[k,,]) <- 0.5 + 0.5*runif(p)
	}
	OmegaINV <- diag(q)
	# compute sufficient stats 
	suf_stat <- compute_sufficient_statistics_given_mu(y = y, 
	 z = z, K = K, x_data = x_data, mu = mu)

	v_r <- numeric(p) #indicates the non-zero values of Lambdas
	for( r in 1:p ){
		v_r[r] <- min(r,q)
	}
	T_INV <- array(data = 0, dim = c(p,p))
	diag(T_INV) <- diag(var(x_data))
	diag(T_INV) <- 1/diag(T_INV)
	ksi <- colMeans(x_data)
	priorConst1 <- array(diag(T_INV)*ksi, dim =c(p,1))
	# now simulate mu and Lambda
	f2 <- compute_A_B_G_D_and_simulate_mu_Lambda_CUU(SigmaINV = SigmaINV, 
                suff_statistics = suf_stat, OmegaINV = OmegaINV, 
                K = K, priorConst1 = priorConst1, T_INV = T_INV, v_r = v_r)
	# f2$mu contains the simulated means
	# f2$Lambdas contains the simulated factor loadings

}


