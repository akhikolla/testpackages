dblpm_mcmc <- function(network, niter, burnin, thin, x_var, w_var, gamma_var, beta_var, verbose = T)
{
  cpp_dblpm_mcmc(network$dim$tframes, network$dim$N, network$dim$M, network$dim$L, network$dim$D,
                 network$edgelist-1,
                 network$like_pars$x, network$like_pars$w, network$like_pars$gamma, network$like_pars$beta,
                 network$tau_pars$tauw, network$tau_pars$tauw0, network$tau_pars$taugamma, network$tau_pars$taugamma0, network$tau_pars$taubeta, network$tau_pars$taubeta0,
                 network$hyper_pars$taux, network$hyper_pars$delta, network$hyper_pars$aw, network$hyper_pars$bw, network$hyper_pars$agamma, network$hyper_pars$bgamma, network$hyper_pars$abeta, network$hyper_pars$bbeta,
                 niter, burnin, thin,
                 x_var, w_var, gamma_var, beta_var,
                 verbose)
}


