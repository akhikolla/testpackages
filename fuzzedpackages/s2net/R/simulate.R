# simulated data

simulate_extra <- function(n_source = 100, n_target = 100, p = 1000, shift = 10, scenario = "same", response = "linear", sigma2 = 2.5){
  xL = matrix(rnorm(n_source*p, 0, 0.4), ncol = p)
  
  switch (scenario,
          same = {
            beta_lucky = 5*c(rep(1, 5),rep(-1, 5), rep(0, p - 10))/sqrt(10)
            
            xU = matrix(rnorm(n_target*p, 0, 0.4), ncol = p)
            etaL = xL %*% beta_lucky
            etaU = xU %*% beta_lucky 
          },
          lucky = {
            beta_lucky = 5*c(rep(1, 5),rep(-1, 5), rep(0, p - 10))/sqrt(10)
            
            xU = matrix(rnorm(n_target*p, 0, 0.4), ncol = p)
            xU[,1:10] = xU[,1:10] + shift
            
            etaL = xL %*% beta_lucky
            etaU = xU %*% beta_lucky
          },
          unlucky = {
            beta_unlucky = 5*c(rep(1, 10), rep(0, p - 10))/sqrt(10)
            
            xU = matrix(rnorm(n_target*p, 0, 0.4), ncol = p)
            xU[,1:10] = xU[,1:10] + shift
            
            etaL = xL %*% beta_unlucky
            etaU = xU %*% beta_unlucky
          }
  )
  
  switch (response,
          linear = {
            yL = etaL +  rnorm(n_source, 0, sqrt(sigma2))
            yU = etaU +  rnorm(n_target, 0, sqrt(sigma2))
          },
          logit = {
            yL = as.factor(rbinom(n_source, 1, 1/(1 + exp(-etaL))))
            yU = as.factor(rbinom(n_target, 1, 1/(1 + exp(-etaU))))
          }
  )
  
  return(
    list(
      xL = xL,
      yL = yL,
      xU = xU,
      yU = yU
    )
  )
}

x_2_grp = function(N, p1, var1, cor1, p2, var2, cor2, mu){
  S = matrix(cor1, nrow = p1, ncol = p1)
  diag(S) = var1
  X1 = mvrnorm(N, mu = mu[1:p1], Sigma = S)
  S = matrix(cor2, nrow = p2, ncol = p2)
  diag(S) = var2
  X2 = mvrnorm(N, mu = mu[(p1+1):(p1+p2)], Sigma = S)
  cbind(X1, X2)
}

simulate_groups = function(n_source = 100, n_target = 100, p = 200, response = "linear"){
  
  xL = x_2_grp(n_source, p/2, 1, 0.8, p/2, 0.05, 0.01, rep(0, p)) #cos(seq(0, 4*pi, length.out = p))
  
  xU = x_2_grp(n_target, p/2, .1, 0.01, p/2, 1, 0.5, rep(0, p)) #sin(seq(0, 4*pi, length.out = p))
  
  beta.S = rep(0, p)
  beta.S[sample(1:floor(p/2), 5)] = 1
  beta.S[sample(ceiling(p/2 + 1):p, 5)] = 1
  
  beta.T = beta.S * runif(p, 0.9, 1.1)
  
  switch (response,
    linear = {
      SNR = 4 # signal-to-noise ratio
      yL = xL %*% beta.S
      k = sd(yL)/sqrt(SNR)
      yL = yL + k*rnorm(n_source) - mean(yL)
      
      yU = xU %*% beta.T
      yU = yU + k*rnorm(n_target) - mean(yL)
    },
    logit = {
      yL = xL %*% beta.S
      yL = as.factor(rbinom(n_source, size = 1, prob =  (1 + exp(-yL))^-1))
      
      yU = xU %*% beta.T
      yU = as.factor(rbinom(n_target, size = 1, prob =  (1 + exp(-yU))^-1))
      
    }
  )
  return(
    list(
      xL = xL,
      yL = yL,
      xU = xU,
      yU = yU
    )
  )
  
}

#prompt(simulate_groups, "man/simulate_groups.Rd")

# prompt(simulate_extra, "man/simulate_extra.Rd")

# shootout_data = function(nL.train = 135, nU.train = 200, nU.valid = 20){
# 
#   data("shootout")
#   ids = unique(shootout$id)
#   
#   idx.L.train = sample(ids, nL.train)
#   ids = setdiff(ids, idx.L.train)
#   idx.U.train = sample(ids, nU.train)
#   ids = setdiff(ids, idx.U.train)
#   
#   idx.U.valid = sample(ids, nU.valid)
#   ids = setdiff(ids, idx.U.valid)
#   
#   mean.yL.train = mean(shootout$label[shootout$id %in% idx.L.train & shootout$instrument=="1"])
#   
#   data = list(
#     train = list(
#       xL = as.matrix(shootout[shootout$id %in% idx.L.train & shootout$instrument=="1", -(1:4)]),
#       yL = as.matrix(shootout$label[shootout$id %in% idx.L.train & shootout$instrument=="1"] - mean.yL.train),
#       xU = as.matrix(shootout[shootout$id %in% idx.U.train & shootout$instrument=="2", -(1:4)])
#     ),
#     valid = list(
#       xU = as.matrix(shootout[shootout$id %in% idx.U.valid & shootout$instrument=="2", -(1:4)]),
#       yU = as.matrix(shootout$label[shootout$id %in% idx.U.valid & shootout$instrument=="2"] - mean.yL.train)
#     ),
#     test = list(
#       xU = as.matrix(shootout[shootout$id %in% c(idx.U.train, ids) & shootout$instrument=="2", -(1:4)]),
#       yU = as.matrix(shootout$label[shootout$id %in% c(idx.U.train, ids) & shootout$instrument=="2"] - mean.yL.train)
#     )
#   )
#   return(data)
# }

# prompt(shootout_data, "man/shootout_data.Rd")
