## the minus log-likelihood function

lik_con_1d <-
function(param, m0, r0, tau, yij, delta, tij, n_j, t0)
{

nll <- .Call('mloglik', param, m0, r0, tau, yij, delta, tij, n_j, t0)

return(nll)
}

## the gradient of the minus log-likelihood function

gr_con_1d <-
function(param, m0, r0, tau, yij, delta, tij, n_j, t0)
{
	d <- .Call('devlik', param, m0, r0, tau, yij, delta, tij, n_j, t0)
	
	return(d)
}
