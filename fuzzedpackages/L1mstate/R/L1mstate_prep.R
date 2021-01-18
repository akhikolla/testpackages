#############################################
#####  Prepare long-format data         #####
#############################################

#############################################
#####  Information needed for Cox model #####
#############################################

coxinfo=function(x, y){
  
  N0=nrow(x)
  oi=order(y[, "status"], decreasing=TRUE)
  x=x[oi, ];y=y[oi, ]
  oi=order(y[, "time"])
  x=x[oi, ];y=y[oi, ]
  
  ## remove the first censored cases
  i1=which(y[, "status"]==1)
  if(length(i1)==0){
    x=x; y=y
  } else {
    mi1=min(i1)-1
    if (mi1!=0) {
      x=x[-c(1:mi1), ];y=y[-c(1:mi1), ]
    }}
  ty=y[, "time"];tevent=y[, "status"]
  N=nrow(x);n1=sum(y[, "status"])
  
  dty=duplicated(ty) # ties
  
  ### for calculation of log-likelihood
  if (any(dty)) {
    tevent0=tevent
    tevent0[which(dty)]=0
    
    ievent=cumsum(tevent0);loc1=which(tevent0==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=tapply(tevent==1, ievent, sum)
  } else {
    ievent=cumsum(tevent);loc1=which(tevent==1)
    nevent=table(ievent);n=length(unique(ievent))
    nevent1=rep(1, n)
  }
  x = as.matrix(x)
  
  return(list(x=x, N0=N0, tevent=tevent, N=N, nevent=nevent, nevent1=nevent1, loc1=loc1, n=n))
}

#############################################
#####  Prepare data for log-likelihood  #####
#############################################

l1mstateprep = function(long_dt){
  
  Q=length(unique(long_dt$trans))
  
  dt = long_dt[long_dt$trans == 1,] 
  drops = c("id","from","to","Tstart","Tstop","time","status","trans")
  x = dt[,!(names(dt) %in% drops)]
  keep = c("time", "status")
  y = dt[,(names(dt) %in% keep)]
  data1 = coxinfo(x,y)
  data = list(data1)
  
  if(Q>1){
    for(q in (2:Q)){
      dt = long_dt[long_dt$trans == q,] 
      drops = c("id","from","to","Tstart","Tstop","time","status","trans")
      x = dt[,!(names(dt) %in% drops)]
      keep = c("time", "status")
      y = dt[,(names(dt) %in% keep)]
      dataq = coxinfo(x,y)
      data[[q]] = dataq
    }
  }
  
  return(data)
}