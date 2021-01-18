
actgt <- function(y, X, xs, hyps, maxit=0, alpha=0.05){# X here is without intercept
  sqrW = sqrt(mean(y)*(1-mean(y)) ) # sqrt of covariance of y
  
  WIHZ = sqrW *(sweep(X,2,colMeans(X))) ## W^{1/2}*(I-H)Z
  IHZ = WIHZ/sqrW  ##(I-H)Z
  if(maxit==0){
    # the full model
    Tf = sum(colSums(y*IHZ[,xs,drop=FALSE])^2)# full test statistic = y^t (I-H) Z Z^t (I-H) y
    Lamf = round(eigen(tcrossprod(WIHZ[,xs,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8) # eigen(W^{1/2}*(I-H)Z%*%Z^t%*%(I-H)%*%W^{1/2}) 
    Cf = criticalvalue(Lamf,alpha=alpha )
    
    # the selected model
    Ts = sum(colSums(y*IHZ[,hyps,drop=FALSE])^2)
    Lams = round(eigen(tcrossprod(WIHZ[,hyps,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8)  
    Cs = criticalvalue(Lams,alpha=alpha  )
    
    ### the approximate procedure
    if(Tf< Cf || Ts< Cs){ 
      rid = "Not Reject"
    } else if(Ts > Cf){
      rid = "Reject"
    } else{
      rvs = setdiff(xs, hyps) 
      it = sapply(rvs,function(i)  sum(colSums(y*IHZ[,i,drop=FALSE])^2)   ) ## individual test ststistics
      iw = sapply(rvs,function(i) sum(WIHZ[,i,drop=FALSE]^2) ) # individual weights 
      tw = sort(it/iw) 
      TW = names(tw)
      riw = iw[TW]
      rit = it[TW]
      Ls = sum(Lams) + c(0, cumsum(riw) )  ## levels for tilde(R)
      Tlmin = c(Ts) + c(0,cumsum(rit) )   ## tmins for tilde(R)
      
      if(!tactrue(Tlmin,hyps,Cf,WIHZ,alp=alpha )){
        rid = "not reject"
      } else if(tacmax(Tlmin,Ls,tw,Cf,Lamf,Lams, alp=alpha )){
        rid = "reject"
      } else{
        rid = "unsure"
      }
      
    }
    if(rid == "Reject" || rid=="reject"){
      rid="reject"
    } else if(rid == "not reject" || rid == "Not Reject" ){
      rid = "not reject"
    } else{
      rid = "unsure"
    }
    rr = c(rid, 0)
    names(rr)=c("Result", "Iterations")
    return(rr)
    
  } else{
    # ig = sapply(xs,function(i)  sum(colSums(y*IHZ[,i,drop=FALSE])^2)   ) ## individual test ststistics
    # ie = sapply(xs,function(i) sum(WIHZ[,i,drop=FALSE]^2) ) # individual weights 
    rr = unlist(actgt_it(y=y,Tmatrix = IHZ,Cmatrix = WIHZ,fxs = xs,sxs=hyps,maxIt = maxit,a = alpha ) )
    names(rr)=c("Result", "Iterations")
    if(rr[1] == "Reject" || rr[1]=="reject"){
      rr[1]="reject"
    } else if(rr[1] == "not reject" || rr[1] == "Not Reject" ){
      rr[1] = "not reject"
    } else{
      rr[1] = "unsure"
    }
    rr[2] = as.numeric(rr[2])-1 ## count the iterations
    return(rr)
    
  }
  
  
}




actgt_it <- function(y, Tmatrix, Cmatrix, fxs, sxs,Tf,Lamf,Cf,Ts,Lams,Cs, count=1, maxIt=1,a=0.05){
  ## the full model
  if(missing(Tf)|| missing(Lamf)|| missing(Cf) ){
    Tf = sum(colSums(y*Tmatrix[,fxs,drop=FALSE])^2) # full test statistic = y^t (I-H) Z Z^t (I-H) y
    Lamf = round(eigen(tcrossprod(Cmatrix[,fxs,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8) #  w^{1/2}*(I-H)Z%*%t(  )
    Cf = criticalvalue(Lamf, alpha = a   )# RP method in c++ 
  }
  
  ## the selected model
  if(missing(Ts)|| missing(Lams)|| missing(Cs)){
    Ts = sum(colSums(y*Tmatrix[,sxs,drop=FALSE])^2)
    Lams = round(eigen(tcrossprod(Cmatrix[,sxs,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8)  
    Cs = criticalvalue(Lams,alpha = a   )
  }
  
  ### the approximate procedure
  if(Tf< Cf || Ts< Cs){
    rid = "Not Reject"
  } else if(Ts > Cf){
    rid = "Reject"
  } else{
    rvs = setdiff(fxs, sxs) # fxs - sxs
    it = sapply(rvs,function(i)  sum(colSums(y*Tmatrix[,i,drop=FALSE])^2)   )  ## individual test statistics
    iw = sapply(rvs,function(i) sum(Cmatrix[,i,drop=FALSE]^2) ) ## individual weights
    tw = sort(it/iw) ## sort the ratios of it to iw 
    TW = names(tw) ## get the names of sorted individuals
    riw = iw[TW]  
    rit = it[TW]  
    Ls = sum(Lams) + c(0, cumsum(riw) )  ## levels for tilde(R)
    Tlmin = c(Ts) + c(0,cumsum(rit) )   ## tmins for tilde(R)
    
    if(!tactrue(Tlmin,sxs,Cf,Cmatrix,alp=a  )){
      rid = "not reject"
    } else if(tacmax(Tlmin,Ls,tw,Cf,Lamf,Lams,alp=a   )){
      rid = "reject"
    } else{  
      rid = "unsure"
    }
  }
  
  if(count>maxIt) rid="UNSURE" 
  
  if(rid=="unsure"){
    RMVs = names(sort(it,decreasing = TRUE))
    res1 = actgt_it(y=y, Tmatrix = Tmatrix, Cmatrix = Cmatrix, fxs = setdiff(fxs, RMVs[1]),sxs,Ts=Ts, Lams = Lams, Cs=Cs,count = count+1,maxIt=maxIt,a=a  )
    rid = res1$res # go with f-, s
    count <- res1$count
    if(rid =="reject" || rid == "Reject"){
      res2 = actgt_it(y=y, Tmatrix = Tmatrix, Cmatrix = Cmatrix, fxs,sxs = c(sxs, RMVs[1]),Tf = Tf, Lamf=Lamf,Cf = Cf,count=count+1,maxIt=maxIt,a=a ) ## go with f, s+
      rid = res2$res 
      count = res2$count
      
    }
  } 
  
  
  return(list(res = rid, count = count))
}



##check if tmin is above cmax 
tacmax = function(tmins,levels, tw, cf,lf,ls,alp ){
  # levels are for which tmin calculated
  # tw = ratios of teststatistic/weights
  # lf: upper bound of eigenvalue vector
  # ls: lower bound of eigenvalue vector
  above =TRUE
  C = cf
  l1 = levels[length(levels)]
  l2=0
  
  while(tmins[1] <= C ){
    ii = max(which(tmins<C))
    l2 = (C-tmins[ii])/tw[ii] + levels[ii]
    C = criticalvalue(getL(lf,ls,l2),alpha = alp)
    if(tmins[1] > C){
      break
    } else if (l1-l2 > 1e-4){
      l1 = l2
    } else{
      above=FALSE
      break
    }
  }
  
  return(above)
}

##check if tmin is above ctrue for ~R
tactrue = function(tmins, hyxs, cfull, Wmatrix,alp ){ 
  # tmins is the minimal test statistics corresponding to each level between submodel and full model
  # hyxs is the selected set
  # cfull is the starting point
  #Wmatrix is used to calculate the critical value
  above=TRUE 
  C = cfull 
  
  while(tmins[1] < C ){
    ii = max(which(tmins<C))
    if(ii>1) chyp = c(hyxs, names(tmins[2:ii]) ) else chyp= hyxs
    C = criticalvalue(round(eigen(tcrossprod(Wmatrix[,chyp,drop=FALSE]),symmetric = TRUE,only.values = TRUE)$values,8) ,alpha = alp  )
    if(tmins[ii] <= C){
      above = FALSE
      break
    }
  }
  return(above)
}


