
SBTMProbs <- function(adj_cube, allocations)
{
  z <- CollapseLabels(allocations)
  tframes <- dim(adj_cube)[3]
  N <- dim(adj_cube)[2]
  K <- max(z)
  K1 <- K + 1
  
  x <- adj_cube
  y <- array(0,c(N,N,tframes))
  for (t in 1:tframes) for (i in 1:N) for (j in 1:N) if (z[t,i] != 0) if (z[t,j] != 0) y[i,j,t] = 1
  
  pi_probs <- pi_counts <- matrix(0,K1,K1)
  pi_total_counts <- rep(0,K1)
  for (t in 2:tframes) for (i in 1:N) 
  {
    g <- z[t-1,i] + 1
    h <- z[t,i] + 1
    pi_counts[g,h] = pi_counts[g,h] + 1
    pi_total_counts[g] = pi_total_counts[g] + 1
  }
  for (g in 1:K1) for (h in 1:K1) pi_probs[g,h] = pi_counts[g,h] / pi_total_counts[g]
  
  eta_counts <- zeta_counts <- u00_counts <- u01_counts <- u10_counts <- u11_counts <- matrix(0,K,K)
  for (i in 1:(N-1)) for (j in (i+1):N) if (y[i,j,1] > 0)
  {
    g <- z[1,i]
    h <- z[1,j]
    if (g > 0) if (h > 0)
    {
      eta_counts[g,h] = eta_counts[g,h] + x[i,j,1]
      if (g != h) eta_counts[h,g] = eta_counts[h,g] + x[i,j,1]
      zeta_counts[g,h] = zeta_counts[g,h] + 1 - x[i,j,1]
      if (g != h) zeta_counts[h,g] = zeta_counts[h,g] + 1 - x[i,j,1]
    }
  }
  for (t in 2:tframes) for (i in 1:(N-1)) for (j in (i+1):N) if (y[i,j,t] > 0)
  {
    g <- z[t,i]
    h <- z[t,j]
    if (g > 0) if (h > 0)
    {
      if (y[i,j,t-1] == 0) 
      {
        eta_counts[g,h] = eta_counts[g,h] + x[i,j,t]
        if (g != h) eta_counts[h,g] = eta_counts[h,g] + x[i,j,t]
        zeta_counts[g,h] = zeta_counts[g,h] + 1 - x[i,j,t]
        if (g != h) zeta_counts[h,g] = zeta_counts[h,g] + 1 - x[i,j,t]
      }
      else
      {
        u00_counts[g,h] = u00_counts[g,h] + (1-x[i,j,t-1])*(1-x[i,j,t])
        u01_counts[g,h] = u01_counts[g,h] + (1-x[i,j,t-1])*(x[i,j,t])
        u10_counts[g,h] = u10_counts[g,h] + (x[i,j,t-1])*(1-x[i,j,t])
        u11_counts[g,h] = u11_counts[g,h] + (x[i,j,t-1])*(x[i,j,t])
        if (g != h)
        {
          u00_counts[h,g] = u00_counts[h,g] + (1-x[i,j,t-1])*(1-x[i,j,t])
          u01_counts[h,g] = u01_counts[h,g] + (1-x[i,j,t-1])*(x[i,j,t])
          u10_counts[h,g] = u10_counts[h,g] + (x[i,j,t-1])*(1-x[i,j,t])
          u11_counts[h,g] = u11_counts[h,g] + (x[i,j,t-1])*(x[i,j,t])
        }
      }
    }
  }
  theta_probs <- p_probs <- q_probs <- matrix(0,K,K)
  for (g in 1:K) for (h in 1:K)
  {
    theta_probs[g,h] = eta_counts[g,h] / (eta_counts[g,h] + zeta_counts[g,h])
    p_probs[g,h] = u01_counts[g,h] / (u01_counts[g,h] + u00_counts[g,h])
    q_probs[g,h] = u10_counts[g,h] / (u10_counts[g,h] + u11_counts[g,h])
  }
  list(Pi = pi_probs, Theta = theta_probs, P = p_probs, Q = q_probs)
}
