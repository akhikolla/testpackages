slpm_nga <- function(X, K, var_pars_init, hyper_pars = NULL, tol = 0.01, n_iter_max = 100000, natural_gradient = T, learning_rate_factor_up = 2, learning_rate_factor_down = 2, verbose = F)
{
  M <- nrow(X)
  N <- ncol(X)
  if (is.null(hyper_pars)) hyper_pars = list(delta = rep(0.001,K), a_gamma = rep(1,K), b_gamma = rep(1,K))
  
  res <- cpp_SLPM_Optimisation(adj = X, 
                        var_alpha_u = var_pars_init$alpha_u_tilde, var_alpha_v = var_pars_init$alpha_v_tilde, 
                        var_beta_u = var_pars_init$beta_u_tilde, var_beta_v = var_pars_init$beta_v_tilde, 
                        var_lambda = var_pars_init$lambda_tilde, 
                        var_delta = var_pars_init$delta_tilde, var_a = var_pars_init$a_tilde, var_b = var_pars_init$b_tilde, 
                        delta = hyper_pars$delta, a_gamma = hyper_pars$a_gamma, b_gamma = hyper_pars$b_gamma, 
                        tol = tol, n_iter_max = n_iter_max, natural_gradient = natural_gradient, 
                        learning_rate_factor_up = learning_rate_factor_up, learning_rate_factor_down = learning_rate_factor_down, verbose = verbose)
  
  list(computing_time = res$computing_time,
       var_pars = list(alpha_u_tilde = res$var_alpha_u, alpha_v_tilde = res$var_alpha_v, 
                       beta_u_tilde = res$var_beta_u, beta_v_tilde = res$var_beta_v, 
                       lambda_tilde = res$var_lambda, 
                       delta_tilde = as.numeric(res$var_delta), a_tilde = as.numeric(res$var_a), b_tilde = as.numeric(res$var_b)),
       learning_rates_u = res$learning_rates_u, learning_rates_v = res$learning_rates_v, 
       elbo_values = as.numeric(res$elbo_values), elbo_init = res$elbo_init, elbo_final = res$elbo_final)
}
