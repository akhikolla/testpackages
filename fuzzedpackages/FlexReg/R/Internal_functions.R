#' @title var.fun
#' @export
#' @keywords internal
#'

var.fun <- function(model, mu.chain, phi.chain){
  model.name <- model$model@model_name
  posterior <- model$model

  if("Beta" %in% substr(model.name,1,4)){
    var1 <- NULL; var2 <- NULL
    if (is.na(dim(phi.chain)[2])){
      variance <- apply(mu.chain*(1-mu.chain),2, function(x) x/(1+phi.chain))
    } else  variance <- (mu.chain*(1-mu.chain))/(1+phi.chain)
  }

  if("VIB" %in% substr(model.name,1,3)){
    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    k.chain <- rstan::extract(posterior, pars="k", permuted=T)[[1]]

    if (is.na(dim(phi.chain)[2])){
      var1 <- apply(mu.chain*(1-mu.chain),2, function(x) p.chain*x/(1+phi.chain*k.chain))
      var2 <- apply(mu.chain*(1-mu.chain),2, function(x) (1-p.chain)*x/(1+phi.chain))
      variance <- var1+var2
    } else {
      var1 <- (mu.chain*(1-mu.chain))/apply(phi.chain, 2, function(x) (1+x*k.chain)/p.chain)
      var2 <- apply((mu.chain*(1-mu.chain))/(1+phi.chain), 2, function(x) (1-p.chain)*x)
      variance <- var1+var2
    }
  }

  if("FB" %in% substr(model.name,1,2)){
    p.chain <- rstan::extract(posterior, pars="p", permuted=T)[[1]]
    w.chain <- rstan::extract(posterior, pars="w", permuted=T)[[1]]
    wtilde.chain <- apply(mu.chain, 2, function(x) w.chain*pmin(x/p.chain, (1-x)/(1-p.chain)))
    lambda1.chain <- mu.chain+ apply(wtilde.chain, 2, function(x) (1-p.chain)*x)
    lambda2.chain <- mu.chain+ apply(wtilde.chain, 2, function(x) p.chain*x)
    if (is.na(dim(phi.chain)[2])){
      var1 <- apply(lambda1.chain*(1-lambda1.chain),2, function(x) p.chain*x/(1+phi.chain))
      var2 <- apply(lambda2.chain*(1-lambda2.chain),2, function(x) (1-p.chain)*x/(1+phi.chain))
      variance <- apply(mu.chain*(1-mu.chain)+ apply(wtilde.chain, 2, function(x) x^2*phi.chain*p.chain*(1-p.chain)), 2, function(x) x/(1+phi.chain))
    } else {
      var1 <- (lambda1.chain*(1-lambda1.chain))/(1+phi.chain)
      var2 <- (lambda2.chain*(1-lambda2.chain))/(1+phi.chain)
      variance <- (mu.chain*(1-mu.chain)+apply(wtilde.chain^2*phi.chain,2, function(x) x*p.chain*(1-p.chain)))/(1+phi.chain)
    }
  }
  return(list(variance=variance, var1=var1, var2=var2))
}

#' @title mu.chain.nd
#' @export
#' @keywords internal
#'
mu.chain.nd <- function(posterior, newdata, link.mu){
  beta.chain <- rstan::extract(posterior, pars="beta", permuted=T)[[1]]
  eta.chain <- beta.chain %*% t(newdata)
  if(link.mu == "logit") mu.chain <- apply(eta.chain,c(1,2), function(x) 1/(1+exp(-x))) else
    if(link.mu == "probit") mu.chain <- apply(eta.chain,c(1,2), function(x) pnorm(x)) else
      if(link.mu == "cloglog") mu.chain <- apply(eta.chain,c(1,2), function(x) 1-exp(-exp(x))) else
        if(link.mu == "loglog") mu.chain <- apply(eta.chain,c(1,2), function(x) exp(-exp(x)))
  return(mu.chain)
}

#' @title phi.chain.nd
#' @export
#' @keywords internal
#'
phi.chain.nd <- function(posterior, newdata, link.phi){
  if(link.phi == "identity") {
    phi.chain <- rstan::extract(posterior, pars="phi", permuted=T)[[1]]} else {
      psi.chain <- rstan::extract(posterior, pars="psi", permuted=T)[[1]]
      eta.chain <- psi.chain %*% t(newdata)
      if(link.phi == "log") phi.chain <- apply(eta.chain,c(1,2), function(x) exp(x)) else
        if(link.phi == "sqrt") phi.chain <- apply(eta.chain,c(1,2), function(x) x^2)
    }
  return(phi.chain)
}
