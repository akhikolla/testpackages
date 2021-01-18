corlatent <- function(data, accuracy, n, R, p, lambda1, lambda2, lambda3, distribution = "Gaussian", rule = "AND"){
  C_plus <- generate_Cplus(R,n)
  d <- matrix(0, nrow = p, ncol = p)
  theta_hat <- matrix(0, nrow = p-1, ncol = p)
  alpha_hat <- matrix(0, nrow = p, ncol = p)
  delta_hat <- matrix(0, nrow = n*R, ncol = p)
  p.fac <- rep(1, n*R)
  p.fac[(n*(R - 1)+1):(n*R)] <- 0

  if (class(data)[1] == "matrix") {
    dat <- data
  }

  if (class(data)[1] == "list") {
    dat <- data[[1]]
    for (i in 2:n) {
      dat <- rbind(dat, data[[i]])
    }
  }

  if (class(data)[1] == "array") {
    dat <- data[,,1]
    for (i in 2:n) {
      dat <- rbind(dat, data[,,i])
    }
  }

  if (class(data)[1] == "data.frame") {
    dat <- as.matrix(data)
  }

  Xj_star <- Xj_star_stacking(dat,n,R,p)


  if (distribution == "Gaussian") {
    for (j in 1:p) {
      eeee <- BCD_Gaussiaon_j(j, dat, Xj_star, C_plus, p.fac, accuracy, lambda1, lambda2, lambda3, theta_hat[,j], alpha_hat[,j], delta_hat[,j], R, n, p)
      theta_hat[,j] <- eeee$theta_hat_j
      alpha_hat[,j] <- eeee$alpha_hat_j
      delta_hat[,j] <- eeee$delta_hat_j
    }
  }

  if (distribution == "Ising") {
    for(j in 1:p){
      eeee <- BCD_Ising_j(j, dat, Xj_star, C_plus, p.fac, accuracy, lambda1, lambda2, lambda3, theta_hat[,j], alpha_hat[,j], delta_hat[,j], R, n, p)
      theta_hat[,j] <- eeee$theta_hat_j
      alpha_hat[,j] <- eeee$alpha_hat_j
      delta_hat[,j] <- eeee$delta_hat_j
    }
  }

  omega <- theta_squarematrix(theta_hat,p,rule)
  
  res <- omega
  # AND rule
  if (rule == "AND") {
    res <- res!=0
    res <- res*1
    res <- res + t(res)
    res[which(res==1)] <- 0
    res <- res/2
  }
  # OR rule
  if (rule == "OR") {
    res <- (res+t(res))/2
  }
  res <- res!=0
  res <- res*1
  return(list(omega=omega, theta = res, penalties = c(lambda1, lambda2, lambda3)))
}
