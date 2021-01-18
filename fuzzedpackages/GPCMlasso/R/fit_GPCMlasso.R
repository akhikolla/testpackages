fit_GPCMlasso <- function(model = model, loglik_fun = loglikPCMlasso, score_fun = scorePCMlasso,
                          log_score_fun = loglikscorePCMlasso,
                          design_list = design_list, 
                          control = control, start = NULL, scale_fac = 1, main.effects = TRUE){

  ## start with-expression to use everything from design_list and control
  with(c(design_list,control),{
  
    ## get nodes and weights for Gauss-Hermite quadrature
    her_poly <- gauss.quad(Q, "hermite")
    GHnodes <- her_poly$nodes
    GHweights <- her_poly$weights * exp(GHnodes^2) * dnorm(GHnodes)

  ## lower bounds for variance parameters
  l.bound <- c(rep(-Inf,px-n_sigma), rep(0,n_sigma))

  if(!is.null(lambda)){
    l.lambda <- NA
  }  

  ## get values for penalty weight, start values and lambda grid    
  help_me <- help_fit(model, Y, l.lambda, start, loglik_fun, score_fun, log_score_fun, adaptive,
                      Q,q,I,n,m,response,design,designX, px,
                      GHweights, GHnodes, acoefs, lambda2, cvalue, n_sigma,
                      l.bound, trace, log.lambda, weight.penalties, scale_fac, ada.lambda,
                      lambda.min, ada.power, cores, null_thresh, DSF, gradtol, iterlim, steptol,
                      main.effects = main.effects)
  
  weight <- help_me$weight
  alpha.start <- help_me$alpha.start
  if(is.null(lambda)){
    lambda <- help_me$lambda
  }
  
  ## start estimation along lambda
  coef.final <- coef.orig <- matrix(0, nrow=length(lambda), ncol = px)
  logLik <- df <- c()
  
  ########################
  for(l in seq_along(lambda)){
    if(trace){
    cat(paste0(l,". lambda out of ", length(lambda),":"), lambda[l],"\n")
    }
suppressWarnings(
    m.opt <- try(nlm(log_score_fun, alpha.start,
                        Q = Q, q = q, I = I, n = n, Y = response, X = design, Z = designX, px = px,
                        GHweights = GHweights, GHnodes = GHnodes,
                        acoefs = acoefs, lambda = lambda[l], scale_fac = scale_fac,
                        lambda2 = lambda2, cvalue = cvalue, cores = cores, weight = weight,  n_sigma = n_sigma,
                         gradtol = gradtol, iterlim = iterlim, check.analyticals = FALSE, 
                        steptol = steptol)))

    coefs.l <- try(m.opt$estimate)

    if(inherits(m.opt, "try-error")){
        m.opt <- try(nlminb(start = alpha.start, objective = loglik_fun, gradient = score_fun,
               Q = Q, q = q, I = I, n = n, Y = response, X = design, Z = designX, px = px,
               GHweights = GHweights, GHnodes = GHnodes,
               acoefs = acoefs, lambda = lambda[l], scale_fac = scale_fac,
               lambda2 = lambda2, cvalue = cvalue, cores = cores, weight = weight,  n_sigma = n_sigma,
               control=list(eval.max=500, iter.max=500, step.min=0.01),
               lower=l.bound))
        coefs.l <- m.opt$par
    }
    alpha.start <- coefs.l
      
    coefs.l[abs(coefs.l)<null_thresh & rowSums(abs(acoefs))!=0] <- 0
      
    coef.orig[l,] <- coefs.l
    
      coefs.l <- round(coefs.l, precision)

      if(sum(colSums(abs(acoefs))!=0)==m*I){
        df.l <- sum(abs(t(acoefs)%*%coefs.l)!=0)
      }else{
        df.l <- 0
        coefs.pen <- coefs.l[rowSums(abs(acoefs))!=0]
        start.coefs <- 1
        for(i in 1:I){
          coefs.i <- coefs.pen[start.coefs:(start.coefs-1+q[i]*m)]
          gamma.index <- matrix(1:(q[i]*m),nrow=m)
          for(ii in 1:m){
            coefs.ii <- coefs.i[gamma.index[ii,]]
            df.l <- df.l + sum(unique(coefs.ii)!=0)
          }
          start.coefs <- q[i]*m+start.coefs
        }

      }
      df.l <- df.l + sum(rowSums(abs(acoefs))==0)
      
      logLik.l <- -loglik_fun(coefs.l, Q = Q, q = q, I = I, n = n,
                             Y = response, X = design, Z = designX, px = px,
                             GHweights = GHweights, GHnodes = GHnodes, 
                             acoefs = acoefs, lambda = 0, lambda2 = 0, cvalue = cvalue, 
                             cores = cores, weight = weight,  n_sigma = n_sigma, scale_fac = 1)

      coef.final[l,] <- coefs.l
      logLik[l] <- logLik.l
      df[l] <- df.l
  }

  coef.rescal <- round(t(t(coef.orig)/sd.vec), precision)
  
  
  ret.list <- list(coefficients = coef.final, coef.rescal = coef.rescal,
                   logLik = logLik, lambda = lambda, df = df)
  
  return(ret.list)
  
  ## end of with-expression
  })

}
