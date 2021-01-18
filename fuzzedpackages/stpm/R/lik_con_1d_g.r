## the minus log-likelihood function

lik_con_1d_g <-
function(param, m0, r0, tau, yij, delta, tij, n_j, geno_a, geno_b, geno_q, geno_f, geno_f1, geno_mu, t0)
{

nll <- .Call('mloglik_g_2', param, m0, r0, tau, yij, delta, tij, n_j, geno_a, geno_b, geno_q, geno_f, geno_f1, geno_mu, t0)

return(nll)
}

## the gradient of the minus log-likelihood function

gr_con_1d_g <-
function(param, m0, r0, tau, yij, delta, tij, n_j, geno_a, geno_b, geno_q, geno_f, geno_f1, geno_mu, t0)
{
	d <- .Call('devlik_g_2', param, m0, r0, tau, yij, delta, tij, n_j, geno_a, geno_b, geno_q, geno_f, geno_f1, geno_mu, t0)

	return(d)
}
