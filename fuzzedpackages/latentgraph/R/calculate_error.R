calculate_error <- function(theta_hat, theta_test, alpha_hat, alpha_test, delta_hat, delta_test){

  theta_error <- sum((theta_hat - theta_test)^2)
  alpha_error <- sum((alpha_hat - alpha_test)^2)
  delta_error <- sum((delta_hat - delta_test)^2)

  error <- max(theta_error,alpha_error,delta_error)
  return(error)
}
