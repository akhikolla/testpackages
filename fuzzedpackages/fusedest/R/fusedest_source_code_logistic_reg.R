######## Iterative scheme #############################

fusedest_logit <- function(X = X, y = y, label_dc = label_dc, H = H, rho = rho, no_lambda = no_lambda, lambda_list = lambda_list,
                                         beta_ini = beta_ini, max_iter = max_iter, tol_err = tol_err, no.cores = no.cores){

  ######## Internal functions #######


  Computeb <- function(p, theta, xi, o_edge_ext, deg_dc, ind_edge_strt) {
    .Call('_fusedest_Computeb', PACKAGE = 'fusedest', p, theta, xi, o_edge_ext, deg_dc, ind_edge_strt)
  }

  IRLSLogisticReg <- function(X, y, a, b, beta_ini, max_iter, tol_err) {
    .Call('_fusedest_IRLSLogisticReg', PACKAGE = 'fusedest', X, y, a, b, beta_ini, max_iter, tol_err)
  }

  UpdateAlphal2Norm <- function(theta01, theta02, tau, p, q_H, lambda) {
    .Call('_fusedest_UpdateAlphal2Norm', PACKAGE = 'fusedest', theta01, theta02, tau, p, q_H, lambda)
  }

  UpdateTheta <- function(theta01, alpha, tau, beta, xi, theta02) {
    .Call('_fusedest_UpdateTheta', PACKAGE = 'fusedest', theta01, alpha, tau, beta, xi, theta02)
  }

  UpdateTau <- function(tau01, theta01, theta02, alpha) {
    .Call('_fusedest_UpdateTau', PACKAGE = 'fusedest', tau01, theta01, theta02, alpha)
  }

  UpdateXi <- function(xi01, beta, theta) {
    .Call('_fusedest_UpdateXi', PACKAGE = 'fusedest', xi01, beta, theta)
  }

  ComputePrimalError <- function(p, q_H, theta01, theta02, theta03, theta04, beta01, beta02, alpha) {
    .Call('_fusedest_ComputePrimalError', PACKAGE = 'fusedest', p, q_H, theta01, theta02, theta03, theta04, beta01, beta02, alpha)
  }

  ComputeDualError <- function(p, q_H, xi01, xi02, tau) {
    .Call('_fusedest_ComputeDualError', PACKAGE = 'fusedest', p, q_H, xi01, xi02, tau)
  }

  ComputeLogitLoss <- function(X, y, beta, ind_strt, n_dc) {
    .Call('_fusedest_ComputeLogitLoss', PACKAGE = 'fusedest', X, y, beta, ind_strt, n_dc)
  }

  Suml2Distance <- function(p, q_H, a, b) {
    .Call('_fusedest_Suml2Distance', PACKAGE = 'fusedest', p, q_H, a, b)
  }


  LogitRegRcppWols.apply02 <- function(i){

    ind_i <- c(ind_strt[i]:ind_end[i])
    IRLSLogisticReg(X = X[ind_i,], y = y[ind_i], a = a[i], b = b[i,], beta_ini = rep(0, p), max_iter = max_iter_beta, tol_err = tol_err_beta)$beta
  }


  UpdateBetaLogisticReg <- function(m, p, X, y, a, r, beta_ini,
                                    theta.ij.list, theta.ji.list, xi.ij.list, xi.ji.list, rho,
                                    o_edge_ext, deg_dc, ind_edge_strt, ind_strt, ind_end,
                                    max_iter_beta, tol_err_beta, no.cores){

    LogitRegRcppWols.apply02 <- function(i){

      ind_i <- c(ind_strt[i]:ind_end[i])
      IRLSLogisticReg(X = X[ind_i,], y = y[ind_i], a = a[i], b = b[i,], beta_ini = beta_ini[i,], max_iter = max_iter_beta, tol_err = tol_err_beta)$beta
    }


    beta <- matrix(0, nrow = m, ncol = p)

    if(r == 1){ ##### This sets initial beta with zero!!!
      beta <- beta_ini ####matrix(rnorm(m*p, 0, 0.01), nrow = m, ncol = p)
    }


    if(r  > 1){

      ###### b.list is a 2*q_H x p matrix #######

      b <- Computeb(p = p, theta = c(theta.ij.list, theta.ji.list), xi = c(xi.ij.list, xi.ji.list),
                    o_edge_ext = as.numeric(o_edge_ext), deg_dc = deg_dc, ind_edge_strt = ind_edge_strt)$b*rho

      beta <- t(parallel::mcmapply(LogitRegRcppWols.apply02, c(1:m), mc.cores = no.cores))
    }

    return(beta)

  }


  ####### Set up key quantities #################

  strt.int <- Sys.time()


  mem.int.strt <- gc()[2,2]

  table_label_dc <- table(label_dc)
  n_dc <- as.numeric(table_label_dc) ### Number of data points in each data center
  m <- length(n_dc)               #### Number of data centers
  p <- dim(X)[2]                  #### Number of coefficients in regression model

  ####### Create indices for data set (n indices) ##########################################

  ind_strt <- 1
  ind_end <- n_dc[1]

  if(m > 1){
    ind_strt <- as.numeric(c(1, cumsum(n_dc[1:(m-1)])+1))
    ind_end <- as.numeric(cumsum(n_dc))
  }

  ###### Create indices for comparison graph (m indices) #################################

  q_H <- dim(H)[1]                #### Number of edges in comparison graph H

  deg_dc <- as.numeric(degree(graph_from_edgelist(H, directed = FALSE))) #### igraph function to calculate node degree in H

  ind_edge_strt <- 1

  if(m > 1){
    ind_edge_strt <- c(1, cumsum(deg_dc[1:(m-1)])+1)
  }

  edge_ext <- c(H[,1],H[,2])
  o_edge_ext <- order(edge_ext, decreasing = FALSE) #### A 2*q_H-dimensional vector

  ##### Sort out the dataset!!! ##################################################################

  o.label_dc <- order(label_dc, decreasing = FALSE)
  label_dc02 <- label_dc[o.label_dc]
  X02 <- X[o.label_dc,]
  y02 <- y[o.label_dc]

  ###### Compute quantities that can be re-used during iteration. ################

  a <- rho*deg_dc        ###### It computes "a" in the first line of the iterative scheme #####
  b <- matrix(0, nrow = m, ncol = p)
  max_iter_beta <- 100
  tol_err_beta <- 10^(-6)

  if(is.null(beta_ini)==TRUE){

      beta_ini <- t(parallel::mcmapply(LogitRegRcppWols.apply02, c(1:m), mc.cores = no.cores))
  }

  if(is.null(lambda_list) == TRUE){

    beta_ini.norm <- sqrt(apply(beta_ini^2, 1, sum))
    max.lambda <- max(beta_ini.norm)
    lambda_list <- seq(0.5*max.lambda, 0.01*max.lambda, length = no_lambda)*rho

  }

  mem.int.end <- gc()[2,2]

  end.int <- Sys.time()


  alg.matrix <- matrix(0, nrow = max_iter*no_lambda, ncol = 6)

  iter.counter <- rep(0, no_lambda)
  int.time <- rep(0, no_lambda)
  iter.time <- rep(0, no_lambda)
  other.time <- rep(0, no_lambda)
  obj.val <- rep(0, no_lambda)
  obj.val.erg <- rep(0, no_lambda)
  mem.int.usage <- rep(0, no_lambda)
  mem.iter.usage <- rep(0, no_lambda)

  alpha_list <- vector(mode = "list", length = no_lambda)
  beta_list <- vector(mode = "list", length = no_lambda)

  ###### Compute initial values ################


  for(v in 1:no_lambda){


    mem.iter.strt <- gc()[2,2]

    lambda <- lambda_list[v]

    beta.i.list <- rep(0, q_H*p)
    beta.j.list <- rep(0, q_H*p)
    alpha.ij.list <- rep(0, q_H*p)
    theta.ij.list <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1) #
    theta.ji.list <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1)  #
    xi.ij.list <- rep(0, q_H*p) # #rnorm(q_H*p, 0, 1)
    xi.ji.list <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1)
    tau.ij.list <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1) #

    theta.ij.list.star <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1) #
    theta.ji.list.star <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1)  #
    xi.ij.list.star <- rep(0, q_H*p) # #rnorm(q_H*p, 0, 1)
    xi.ji.list.star <- rep(0, q_H*p) #rnorm(q_H*p, 0, 1)
    tau.ij.list.star <- rep(0, q_H*p)

    beta.erg <- matrix(0, nrow = m, ncol = p)
    alpha.ij.list.erg <- rep(0, q_H*p)

    ##### Run the iterative scheme. Multicores have no advantage when m is small!!!!

    primal.err <- 10^10
    dual.err <- 10^10
    r <- 1


    while(r <= max_iter & primal.err+dual.err >= tol_err){

      strt <- Sys.time()

      ###### Update beta.i #############

      beta <- UpdateBetaLogisticReg(m = m, p = p, X = X, y = y, a = a, r = r, beta_ini = beta_ini,
                                    theta.ij.list = theta.ij.list, theta.ji.list = theta.ji.list,
                                    xi.ij.list = xi.ij.list, xi.ji.list = xi.ji.list, rho = rho,
                                    o_edge_ext = o_edge_ext, deg_dc = deg_dc, ind_edge_strt = ind_edge_strt,
                                    ind_strt = ind_strt, ind_end = ind_end,
                                    max_iter_beta = max_iter_beta, tol_err_beta = tol_err_beta, no.cores = no.cores)

      beta.i.list <- as.vector(t(beta[H[,1],]))
      beta.j.list <- as.vector(t(beta[H[,2],]))

      ###### Update alpha.ij ###########

      alpha.ij.list <- UpdateAlphal2Norm(theta01 = theta.ij.list, theta02 = theta.ji.list,
                                         tau = tau.ij.list, p = p, q_H = q_H, lambda = lambda/rho)$alpha

      ############# Update theta.ij, theta.ji, tau.ij, xi.ij and xi.ji ####

      theta.ij.list.star <- UpdateTheta(theta01 = theta.ji.list, alpha = alpha.ij.list, tau = tau.ij.list, beta = beta.i.list,
                                        xi = xi.ij.list, theta02 = theta.ij.list)
      theta.ji.list.star <- UpdateTheta(theta01 = theta.ij.list, alpha = -alpha.ij.list, tau = -tau.ij.list, beta = beta.j.list,
                                        xi = xi.ji.list, theta02 = theta.ji.list)
      tau.ij.list.star <- UpdateTau(tau01 = tau.ij.list, theta01 = theta.ij.list.star,
                                    theta02 = theta.ji.list.star, alpha = alpha.ij.list)
      xi.ij.list.star <- UpdateXi(xi01 = xi.ij.list, beta = beta.i.list, theta = theta.ij.list.star)
      xi.ji.list.star <- UpdateXi(xi01 = xi.ji.list, beta = beta.j.list, theta = theta.ji.list.star)

      ############ Compute stopping criteria ######

      primal.err <- ComputePrimalError(p = p, q_H = q_H, theta01 = theta.ij.list, theta02 = theta.ji.list,
                                       theta03 = theta.ij.list.star, theta04 = theta.ji.list.star,
                                       beta01 = beta.i.list, beta02 = beta.j.list, alpha = alpha.ij.list)

      dual.err <- ComputeDualError(p = p, q_H = q_H, xi01 = xi.ij.list.star, xi02 = xi.ji.list.star, tau = tau.ij.list.star)

      end <- Sys.time()

      alg.matrix[((v-1)*max_iter + r), 1] <- difftime(end, strt, units="secs")

      alg.matrix[((v-1)*max_iter + r), 3] <- primal.err
      alg.matrix[((v-1)*max_iter + r), 4] <- dual.err



      ########### Compute objective value #########

      l1 <- ComputeLogitLoss(X = X02, y = y02, beta = beta, ind_strt = ind_strt, n_dc = n_dc)
      l2 <-  (lambda/rho)*Suml2Distance(p = p, q_H = q_H, a = alpha.ij.list, b = rep(0, p*q_H))$sum_l2_dist
      obj.val.r <- l1 + l2

      beta.erg <- ((r-1)/r)*beta.erg + (1/r)*beta
      alpha.ij.list.erg <- ((r-1)/r)*alpha.ij.list.erg + (1/r)*alpha.ij.list

      l1.erg <- ComputeLogitLoss(X = X02, y = y02, beta = beta.erg, ind_strt = ind_strt, n_dc = n_dc)
      l2.erg <-  (lambda/rho)*Suml2Distance(p = p, q_H = q_H, a = alpha.ij.list.erg, b = rep(0, p*q_H))$sum_l2_dist
      obj.val.erg.r <- l1.erg + l2.erg

      alg.matrix[((v-1)*max_iter + r), 5] <- obj.val.r
      alg.matrix[((v-1)*max_iter + r), 6] <- obj.val.erg.r

      end.others <- Sys.time()

      alg.matrix[((v-1)*max_iter + r), 2] <- difftime(end.others, strt, units="secs") - difftime(end, strt, units = "secs")


      ########### Replace theta.ij.list with theta.ij.list.star #####

      theta.ij.list <- theta.ij.list.star
      theta.ji.list <- theta.ji.list.star
      tau.ij.list <- tau.ij.list.star
      xi.ij.list <- xi.ij.list.star
      xi.ji.list <- xi.ji.list.star

      #print(c(r, primal.err, dual.err, obj.val))
      r <- r+1
    }

    mem.iter.end <- gc()[2,2]

    iter.counter[v] <- r-1
    int.time[v] <- difftime(end.int, strt.int, units = "secs")
    iter.time[v] <- sum(alg.matrix[c(((v-1)*max_iter + 1):(v*max_iter)),1])
    other.time[v] <- sum(alg.matrix[c(((v-1)*max_iter + 1):(v*max_iter)),2])
    obj.val[v] <- obj.val.r
    obj.val.erg[v] <- obj.val.erg.r

    mem.int.usage[v] <- mem.int.end - mem.int.strt
    mem.iter.usage[v] <- mem.iter.end - mem.iter.strt
    if(v > 1){
      mem.iter.usage[v] <- mem.iter.usage[1]
    }

    alpha_list[[v]] <- alpha.ij.list
    beta_list[[v]] <- beta
  }
  result.list <- list(alg.matrix, iter.counter, int.time, iter.time, other.time, obj.val, obj.val.erg, mem.int.usage, mem.iter.usage, alpha_list, beta_list)
  names(result.list) <- c("alg.matrix", "iter.counter", "int.time", "iter.time", "other.time", "obj.val", "obj.val.erg", "mem.int.usage", "mem.iter.usage", "alpha_list", "beta_list")
  return(result.list)

}
