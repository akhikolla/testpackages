################################
#####  Multi-state model   #####
#####  Lasso (L1-norm)     #####
#####  Coordinate Descent  #####
## Penalized Cross-Validation###
################################

l1mstateR = function(longdt, lambda=NULL, nlambda=100, rlambda=NULL, thresh=1e-7, maxit=1e+5){
  
  N0 = nrow(longdt)
  P = ncol(longdt[,-c(1:8)])
  
  #### process data
  listdt = l1mstateprep(longdt)
  Q = length(listdt)
  
  #### centered and standardized
  xscale = list()
  for(q in (1:Q)){
    xq = listdt[[q]]; x = xq[[1]]
    tem = scaleQ(x)
    xq[[1]] = tem$x; listdt[[q]] = xq
    xscale[[q]] = tem$sd
    rm(tem)
  }
  
  #### lambda path
  if(is.null(lambda)){
    lambda_max = max_lambdaQ(listdt)
    if(is.null(rlambda)){
      rlambda = ifelse(N0>P, 0.0001, 0.01)
    }
    lambda_min = lambda_max*rlambda
    lambda = lambda_max*(lambda_min/lambda_max)^(c(0:(nlambda-1))/(nlambda-1))
  }else{
    nlambda = length(lambda)
  }
  
  #### Main function
  out = l1msQ(listdt, lambda, nlambda, P, thresh, maxit)
  
  nlambdai = out$nlambda
  if(nlambdai==0) return(NULL)
  lambdai = lambda[1:nlambdai]
  
  ## Get the beta without standardizing
  out_betas = out$Beta[1:nlambdai]
  for(i in (1:nlambdai)){
    beta0_i = out_betas[[i]]
    for(q in (1:Q)){
      xscale_q = xscale[[q]]
      for(p in (1:P)){
        beta0_i[q,p] = beta0_i[q,p]/xscale_q[p] 
      }
    }
    out_betas[[i]] = beta0_i
  }
  out$Beta_O = out_betas
  
  temCV = data.frame(lambda=lambdai, cvm=out$ll)
  
  return(list(aBetaSTD=out$Beta, aBetaO=out$Beta_O, 
              fit = temCV, numcovs = P, numtrans = Q))
}
