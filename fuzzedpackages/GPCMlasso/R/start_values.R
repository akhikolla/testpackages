start_values <- function(RSM, GPCM, Y, q, p_pen,p, n_sigma, person_list){
  Y2 <- matrix(NA, ncol=ncol(Y), nrow=nrow(Y))
  for(o in 1:length(person_list)){
    Y2[o,person_list[[o]]+1] <- Y[o,person_list[[o]]+1]
  }
  
  if(!RSM ){
    if(q>1){
      if(GPCM){
        m.ltm <- gpcm(Y2, constraint = "gpcm")
        coef.ltm <- coef(m.ltm)
        sigma.start <- coef.ltm[,q+1]
        coef.ltm <- coef.ltm*sign(sigma.start)
      }else{
        m.ltm <- gpcm(Y2, constraint = "1PL")
        coef.ltm <- coef(m.ltm)
        sigma.start <- coef.ltm[1,q+1]
        coef.ltm <- coef.ltm*sign(sigma.start)
      }
    }else{
      if(GPCM){
        m.ltm <- gpcm(Y2, constraint = "gpcm")
        coef.ltm <- -coef(m.ltm)
        sigma.start <- coef.ltm[,q+1]
        coef.ltm <- coef.ltm*sign(sigma.start)
      }else{
        m.ltm <- rasch(Y2)
        coef.ltm <- -coef(m.ltm)
        sigma.start <- coef.ltm[1,q+1]
        coef.ltm <- coef.ltm*sign(sigma.start)
      } 
    }
    delta.start <- c(t(coef.ltm[,-(q+1)]))
    alpha.start <- c(delta.start,rep(0,p_pen),abs(sigma.start))
  }else{
    alpha.start <- c(rep(0, p-n_sigma), rep(0.1,n_sigma))
    # try(optim(par = alpha.start, fn = loglik_fun, gr = score_fun,
    #           Q = Q, q = q, I = I, n = n, Y = response, X = design, Z = designX, p = p,
    #           GHprobs = GHprobs, GHweights = GHweights, GHnodes = GHnodes,
    #           acoefs = acoefs, lambda = lambda[l],
    #           lambda2 = lambda2, con = con, weight = weight,  n_sigma = n_sigma,
    #           pers_list = person_list, scalefac = scalefac,
    #           lower = l.bound, method = "L-BFGS-B",
    #           control = list(parscale = parscale_now, factr = factr)))
    # 
  }
  
 return(alpha.start)
}