slpm_elbo <- function(X, var_pars, hyper_pars, verbose = F)
{
  cpp_SLPM_ELBO(adj = X, 
                var_alpha_u = var_pars$alpha_u_tilde, var_alpha_v = var_pars$alpha_v_tilde, 
                var_beta_u = var_pars$beta_u_tilde, var_beta_v = var_pars$beta_v_tilde, 
                var_lambda = var_pars$lambda_tilde, 
                var_delta = var_pars$delta_tilde, var_a = var_pars$a_tilde, var_b = var_pars$b_tilde, 
                delta = hyper_pars$delta, a_gamma = hyper_pars$a_gamma, b_gamma = hyper_pars$b_gamma, 
                verbose = verbose)
}

