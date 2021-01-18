##############################################################
###compute the transition intensity matrix Q - a hazard matrix
###based on 'agsurv' function of 'survival' package

cumhaz.l1mstate = function(object, longdt, newdata, cv.method = c("pcvl", "min")){
  ##beta should not STD one
  if(cv.method=="pcvl"){
    beta = object$pBetaO
  }else{
    beta = object$mBetaO
  }
  
  #Compute the baseline hazard and the cumulative baseline hazard 
  bshaz=list()
  varbshaz=list()
  cumbshaz=list()
  haz=list()
  timeQ = list()
  
  for (q in sort(unique(longdt$trans))){
    dt = longdt[longdt$trans == q,] 
    drops = c("id","from","to","Tstart","Tstop","time","status","trans")
    x = dt[,!(names(dt) %in% drops)]
    keep = c("Tstart","Tstop", "status")
    y = dt[,(names(dt) %in% keep)]
    
    #no standardizing x
    x = as.matrix(x)
    nvar = ncol(x)
    
    #### Process
    start = y$Tstart
    stop = y$Tstop
    status = y$status
    wt = rep(1, length(stop))
    time = sort(unique(stop))
    timeQ[[q]] = time
    death = (status==1)
    
    #nevent = status[status==1]
    nevent = as.vector(rowsum(wt*death, stop))
    if(length(nevent)!=0){
      rcumsum = function(x) rev(cumsum(rev(x)))
      score = c(exp(x%*%beta[q,]))
      wrisk = wt*score
      #nrisk = rcumsum(score)
      nrisk = rcumsum(rowsum(wrisk, stop))
      
      #process time
      delta = min(diff(time))/2
      etime = c(sort(unique(start)), max(start)+delta) #unique entry times
      indx = approx(etime, 1:length(etime), time, method = 'constant', rule=2, f=1)$y
      
      #esum = rcumsum(rowsum(score, start)) #not yet entered
      esum = rcumsum(rowsum(wrisk, start))
      nrisk = nrisk - c(esum, 0)[indx]
      
      bshaz[[q]] = nevent/nrisk
      varbshaz[[q]] = nevent/nrisk^2
      cumbshaz[[q]] = cumsum(nevent/nrisk)
      
      haz[[q]] = outer(bshaz[[q]], score, "*")## rows are bshaz - failure time points, columns are risk - # individuals
      ##from here can compute the average hazard value for each time point
      haz[[q]] = data.frame("time"=time, "haz"=rowMeans(haz[[q]]), "trans"=rep(q, length(time)))
    }else{
      haz[[q]] = data.frame("time"=time, "haz"=rep(0, length(time)), "trans"=rep(q, length(time))) 
    }
  }
  
  #Compute the hazard 
  if(missing(newdata)){
    #if newdata missing, return the hazard for original data
    Haz = setNames(do.call(rbind.data.frame, haz), c("time", "Haz", "trans"))
    res = list("time"=timeQ, "baseHaz"=bshaz, "var"=varbshaz, "cumBaseHaz"=cumbshaz, "Haz"=Haz)
    return(res)
  }else{
    haz=list()
    #check whether or not newdata includes trans indicator variable
    if(is.null(newdata$trans)){
      stop("Newdata doesn't have trans variable")
    }else{
      #check whether or not newdata includes the same variables as ordinary data
      ###NEWDATA is not STANDARDIZED
      xnew = newdata[,-1]
      if(ncol(xnew)!=nvar){
        stop("Wrong number of variables in newdata")
      }else{
        for(q in sort(unique(newdata$trans))){
          xnewq = xnew[newdata$trans==q,]
          #### centered and standardized
          xnewq = as.matrix(xnewq)
          
          newrisk = c(exp(xnewq%*%beta[q,]))
          if(length(newrisk)>1){
            haz[[q]] = outer(bshaz[[q]], newrisk, "*")## rows are bshaz - failure time points, columns are risk 
            # individuals
            ##from here can compute the average hazard value for each time point
            haz[[q]] = data.frame("time"=timeQ[[q]], "haz"=rowMeans(haz[[q]]), "trans"=rep(q, length(timeQ[[q]])))
          }else{
            haz[[q]] = bshaz[[q]]*newrisk
            # individuals
            ##from here can compute the average hazard value for each time point
            if(length(haz[[q]])==0){
              haz[[q]] = data.frame("time"=timeQ[[q]], "haz"=rep(0, length(timeQ[[q]])), "trans"=rep(q, length(timeQ[[q]])))
            }else{
              haz[[q]] = data.frame("time"=timeQ[[q]], "haz"=haz[[q]], "trans"=rep(q, length(timeQ[[q]])))
            }
            
          }
        }
        Haz = setNames(do.call(rbind.data.frame, haz), c("time", "Haz", "trans"))
        res = list("time"=timeQ, "baseHaz"=bshaz, "var"=varbshaz, "cumBaseHaz"=cumbshaz, "Haz"=Haz)
        return(res)
      }
    }}
}