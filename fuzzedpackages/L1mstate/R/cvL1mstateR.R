################################
#####  Multi-state model   #####
#####  Lasso (L1-norm)     #####
#####  Coordinate Descent  #####
## Penalized Cross-Validation###
################################

cv.l1mstateR = function(longdt, lambda=NULL, nlambda=100, rlambda=NULL, nfolds=1, foldid=NULL, thresh=1e-7, maxit=1e+5){
  
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
  
  ## Check whether or not to do cross-validation
  if(nfolds==1 & is.null(foldid)){
    fit = data.frame(lambda=lambdai)
    return(list(Beta=out$Beta, BetaO=out$Beta_O, fit=fit))
  }else{
    if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")
    #### Split data for cross-validation
    n = unique(longdt$id)  # get unique id's
    N=length(n)
    foldid=sample(rep(seq(nfolds), length=N))
    
    traintest = list()
    prepk = list()
    for(i in 1:nfolds){
      ni = n[which(foldid!=i)]
      traintest[[i]]=longdt[longdt$id == ni[1],]
      for(j in 2:length(ni)){
        traintest[[i]] = rbind(traintest[[i]],longdt[longdt$id == ni[j],])
      }
      longdt_k=traintest[[i]]
      #### process data
      listdt_k = l1mstateprep(longdt_k)
      Q_k = length(listdt_k)
      #### centered and standardized
      for(q in (1:Q_k)){
        xq = listdt_k[[q]]; x = xq[[1]]
        tem = scaleQ(x)
        xq[[1]] = tem$x; listdt_k[[q]] = xq
        rm(tem)
      }
      prepk[[i]] = listdt_k
    }
    
    #### Cross-validation
    outi = list()
    cvPL = matrix(NA, nrow=nfolds, ncol=nlambdai)
    for(i in 1:nfolds){
      outi[[i]] = cvl1msQ(prepk[[i]], lambdai, nlambdai, P, thresh=1e-7, maxit=1e+5, listdt)
      cvPL[i, 1:outi[[i]]$nlambda] = outi[[i]]$lf[1:outi[[i]]$nlambda] - outi[[i]]$ll[1:outi[[i]]$nlambda]
    }
    
    #### Process results
    cvPL = matrix(cvPL[,1:nlambdai], ncol = nlambdai)
    cvraw = cvPL; nfoldi = apply(!is.na(cvraw), 2, sum); rm(cvPL)
    cvm = apply(cvraw, 2, mean, na.rm = TRUE)
    cvse = sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, mean, na.rm = TRUE)/(nfoldi-1))
    
    indexi = which.max(cvm)
    
    ################################
    ######penalized extensions
    l = rev(lambdai)
    cvl = rev(cvm*N0)
    #lambda maximizes cvl and its corresponding cvl, number of non-zero coefs
    l.cvl = lambdai[indexi]
    cvlMax = cvl[which(l==l.cvl)]
    #degree of freedom = number of non-zeros
    df = rep(0,nlambdai)
    for(j in (1:nlambdai)){
      df[j] = sum(apply(out$Beta[[j]] != 0, 2, sum))  
    }
    nz = rev(df)
    nzMax = nz[which(l==l.cvl)]
    
    nzmin = min(nz)
    cvl0 = max(cvl[which(nz==nzmin)])
    
    if(nzMax != 0){
      #lasso-pcvl: penalized cross-validated log-likelihood
      pcvl = cvl - (((cvlMax - cvl0)/nzMax) *nz)
      l.pcvl = l[which.max(pcvl)]
    }else{
      l.pcvl = l.cvl
    }
    
    indexm = which(lambdai==l.pcvl)
    CV.max = cvm[indexm]
    
    #################################
    ###### result
    temi = rep("", nlambdai)
    temi[indexi] = "cvmax"
    if(indexm==indexi){
      temi[indexm] = "pcvl=cvmax"
    }else{
      temi[indexm] = "pcvl"
    }
    
    temCV = data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, index=temi, stringsAsFactors=FALSE)
    
    return(list(aBetaSTD=out$Beta, aBetaO=out$Beta_O, pBetaSTD=out$Beta[[indexm]], pBetaO=out$Beta_O[[indexm]], mBetaSTD=out$Beta[[indexi]], 
                mBetaO=out$Beta_O[[indexi]], fit=temCV, cvmax=CV.max, lambda.max=lambdai[indexi], lambda.pcvl=lambdai[indexm], 
                numcovs = P, numtrans = Q))
  }
}
