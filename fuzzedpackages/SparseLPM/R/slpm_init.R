slpm_init <- function(X, K, method = "random", threshold = 0.1, stdev = NULL)
{
  X = X + threshold
  M <- nrow(X)
  N <- ncol(X)
  if (is.null(stdev)) stdev <- sqrt(mean(1/X))
  var_pars_init <- list()
  var_pars_init$alpha_u_tilde = matrix(rnorm(M * K, 0, stdev), M, K)
  var_pars_init$alpha_v_tilde = matrix(rnorm(N * K, 0, stdev), N, K)
  var_pars_init$beta_u_tilde = matrix(10 * stdev^2, M, K)
  var_pars_init$beta_v_tilde = matrix(10 * stdev^2, N, K)
  var_pars_init$lambda_tilde = array(1/K, c(M, N, K))
  var_pars_init$delta_tilde = var_pars_init$a_tilde = var_pars_init$b_tilde = rep(1, K)
  if (method == "distance") {
    D_uu <- N^2/(sqrt(X %*% t(X)))
    D_vv <- M^2/(sqrt(t(X) %*% X))
    D_all_1 <- rbind(D_uu,t(1/X))
    D_all_2 <- rbind(D_vv,1/X)
    D_all <- cbind(D_all_1,D_all_2)
    diag(D_all) = 0
    positions <- monoMDS(dist = 0.5*(D_all+t(D_all)), k = K)$points
    var_pars_init$alpha_u_tilde <- positions[1:M,]
    var_pars_init$alpha_v_tilde <- positions[(M+1):(M+N),]
    var_pars_init$beta_u_tilde = matrix(20 * var(as.numeric(var_pars_init$alpha_u_tilde)), M, K)
    var_pars_init$beta_v_tilde = matrix(20 * var(as.numeric(var_pars_init$alpha_v_tilde)), N, K)
  }
  var_pars_init
}
