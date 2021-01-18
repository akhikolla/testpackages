BCD_Ising <- function(dat,accuracy,lambda1,lambda2,lambda3,theta_hat,alpha_hat,delta_hat,R,n,p,Xj_star,C_plus,p.fac){
  #print(lambda1)
  for(j in 1:p){
    eeee <- BCD_Ising_j(j,dat, Xj_star, C_plus, p.fac, accuracy, lambda1, lambda2, lambda3,theta_hat[,j],alpha_hat[,j],delta_hat[,j],R,n,p)
    theta_hat[,j] <- eeee$theta_hat_j
    alpha_hat[,j] <- eeee$alpha_hat_j
    delta_hat[,j] <- eeee$delta_hat_j
    #print(j)
  }
  return(list(theta_hat=theta_hat, alpha_hat=alpha_hat, delta_hat=delta_hat))
}
