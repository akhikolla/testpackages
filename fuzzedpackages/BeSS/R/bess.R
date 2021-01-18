bess = function(x, y, family = c("gaussian", "binomial", "cox"),
                method = "gsection", s.min = 1,
                s.max,
                s.list,
                K.max = 20,
                max.steps = 15,
                glm.max = 1e6,
                cox.max = 20,
                factor = NULL,
                epsilon = 1e-4,
                weights=rep(1,nrow(x)))
{
  family <- match.arg(family)
  if(ncol(x)==1|is.vector(x)) stop("x should be two columns at least!")
  if(missing(family)) stop("Please input family!")
  if(!is.null(factor)) method = "sequential"
  if(family=="binomial")
  {
    if(is.factor(y)){
      y = as.character(y)
    }
    if(length(unique(y))!=2)  stop("Please input binary variable!")else
      if(setequal(y_names<-unique(y),c(0,1))==FALSE)
      {
        y[which(y==unique(y)[1])]=0
        y[which(y==unique(y)[2])]=1
        y=as.numeric(y)
      }
  }
  if(family=="cox")
  {
    if(!is.matrix(y)) y=as.matrix(y)
    if(ncol(y)!=2) stop("Please input y with two columns!")
  }
  if(is.vector(y))
  {
    if(nrow(x)!=length(y)) stop("Rows of x must be the same as length of y!")
  }else{
    if(nrow(x)!=nrow(y)) stop("Rows of x must be the same as rows of y!")
  }


  if(missing(s.max)) s.max=min(ncol(x),round(nrow(x)/log(nrow(x))))

  weights = weights/mean(weights)
  beta0=rep(0,ncol(x))
  if(!is.null(factor)){
    if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x),"g")
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop=FALSE], 2, function(x){
      x = as.factor(x)
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop=FALSE], 2,function(x) {length(unique(x))})-1
    Gi = rep(1:ncol(x), times = group)
    beta0 = rep(beta0, times = group)
    x = model.matrix(~., data = x)[,-1]
  }
  #normalize x
  x=as.matrix(x)
  xs=x
  nm = dim(x)
  n = nm[1]
  p = nm[2]
  one = rep(1, n)
  vn = dimnames(x)[[2]]
  meanx = drop(weights %*% x)/n
  x = scale(x, meanx, FALSE)
  normx = sqrt(drop(weights %*% (x^2)))
  nosignal = normx/sqrt(n) < .Machine$double.eps
  if (any(nosignal))  normx[nosignal] = (.Machine$double.eps) * sqrt(n)
  names(normx) = vn
  x = sqrt(n)*scale(x, FALSE, normx)


  if(method=="sequential"&missing(s.list)) s.list=1:min(ncol(x),round(nrow(x)/log(nrow(x))))

  if(family=="gaussian")
  {
    ys = y
    mu = mean(y*weights)
    y = drop(y - mu)

    #initial beta
    gc()

    if(method=="gsection")
    {
      k = 1
      sL=s.min
      sR=s.max
      beta0R=beta0
      beta0L=beta0

      fit_L1=bess.lm(x,y,
                     beta0=beta0L,
                     s=sL,
                     weights=weights,
                     max.steps=max.steps)

      nullmse=fit_L1$nullmse
      fit_L=fit_L1
      fit_R=bess.lm(x,y,
                    beta0=beta0R,
                    s=sR,
                    weights=weights,
                    max.steps=max.steps)

      beta.fit=cbind(fit_L$beta,fit_R$beta)
      mse=c(fit_L$mse,fit_R$mse)
      lambda=c(fit_L$lambda,fit_R$lambda)
      aic=c(fit_L$AIC,fit_R$AIC)
      bic=c(fit_L$BIC,fit_R$BIC)
      ebic=c(fit_L$EBIC,fit_R$EBIC)

      beta0M=fit_R$beta
      s.list=c(sL,sR)

      while(k<=K.max)
      {
        sM <- round(sL + (sR-sL)*0.618)
        s.list=c(s.list,sM)
        fit_M=bess.lm(x,y,
                      beta0=beta0M,
                      s=sM,
                      weights=weights,
                      max.steps=max.steps)
        cat(k,"-th iteration s.left:",sL," s.split:",sM," s.right:",sR,"\n",sep="")

        beta0M=fit_M$beta
        beta.fit=cbind(beta.fit,beta0M)

        mse=c(mse, fit_M$mse)
        lambda=c(lambda, fit_M$lambda)
        aic=c(aic, fit_M$AIC)
        bic=c(bic, fit_M$BIC)
        ebic=c(ebic, fit_M$EBIC)
        if(abs(fit_L$mse-fit_M$mse)/abs(nullmse*(sM-sL)) > epsilon &
           abs(fit_R$mse-fit_M$mse)/abs(nullmse*(sM-sR)) < epsilon)
        {
          sR <- sM
          fit_R=fit_M
        }else if(abs(fit_L$mse-fit_M$mse)/abs(nullmse) > epsilon &
                 abs(fit_R$mse-fit_M$mse)/abs(nullmse) > epsilon)
        {
          sL <- sM
          fit_L=fit_M
        }else
        {
          sR=sM
          fit_R=fit_M
          sL=s.min
          fit_L=fit_L1
        }

        if(sR-sL==1) break
        fit_ML=bess.lm(x,y,
                       beta0=beta0M,
                       s=sM,
                       weights=weights,
                       max.steps=max.steps)

        fit_MR=bess.lm(x,y,
                       beta0=beta0M,
                       s=sM,
                       weights=weights,
                       max.steps=max.steps)
        if(abs(fit_ML$mse-fit_M$mse)/abs(fit_M$mse) > epsilon &
           abs(fit_MR$mse-fit_M$mse)/abs(fit_M$mse) < epsilon)
        {break}

        k=k+1
      }
    }
    if(method=="sequential")
    {
      nullmse=sum(weights*y^2)/n
      #cat(nullmse,"\\n")
      mse=vector()
      lambda=vector()
      aic=vector()
      bic=vector()
      ebic=vector()
      for(k in 1:length(s.list))
      {
        #cat("select",s.list[k],"variables","\\n")
        if(is.null(factor)){
          if(k == 1){
            fit=bess.lm(x,y,s=s.list[k],weights=weights,
                        max.steps=max.steps,beta0=beta0)
            beta.fit = matrix(fit$beta)
          }else{
            fit=bess.lm(x,y,s=s.list[k],weights=weights,
                        max.steps=max.steps,beta0=beta.fit[,k-1,drop=TRUE])
            beta.fit=cbind(beta.fit,fit$beta)
          }
        }else{
          if(k == 1){
            fit=gbess.lm(x,y,s=s.list[k],weights=weights,Gi=Gi,
                         max.steps=max.steps,beta0=beta0)
            beta.fit = matrix(fit$beta)
            s.list[k] = fit$gr_size
          }else{
            fit=gbess.lm(x,y,s=s.list[k],weights=weights,Gi=Gi,
                         max.steps=max.steps,beta0=beta.fit[,k-1,drop=TRUE])
            beta.fit=cbind(beta.fit,fit$beta)
            s.list[k] = fit$gr_size
          }
        }
        #mse
        mse[k]=fit$mse
        #lambda
        lambda[k]=fit$lambda
        aic[k]=fit$AIC
        bic[k]=fit$BIC
        ebic[k]=fit$EBIC
      }
    }

    beta.fit=sqrt(n)*(beta.fit)/normx
    colnames(beta.fit) = s.list
    rownames(beta.fit) = vn
    coef0=mu-drop(t(beta.fit)%*%meanx)
    names(coef0)=s.list

    xbest=xs[,which(beta.fit[,ncol(beta.fit),drop=TRUE]!=0)]
    bestmodel=lm(ys~xbest, weights=weights)

    out=list(family="bess_gaussian",method=method,beta=beta.fit,coef0=coef0,
             s.list=s.list,meanx=meanx,normx=normx,meany=mu,nsample=n,bestmodel=bestmodel,
             mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic,lambda=lambda,max.steps=max.steps,
             factor=factor)
    class(out)="bess"
    return(out)
  }

  if(family=="binomial")
  {
    beta0=rep(0,ncol(x))
    intercept=0
    gc()

    if(method=="gsection")
    {
      k = 1
      sL=s.min
      sR=s.max
      beta0R=beta0
      beta0L=beta0
      coef0L=intercept
      coef0R=intercept

      fit_L1=bess.glm(x=x,y=y,
                      beta0=beta0L,
                      intercept=coef0L,
                      s=sL,
                      glm.max=glm.max,
                      max.steps=max.steps,
                      weights=weights)

      nulldev=fit_L1$nulldeviance

      fit_L=fit_L1
      fit_R=bess.glm( x=x,y=y,
                      beta0=beta0R,
                      intercept=coef0R,
                      s=sR,
                      glm.max=glm.max,
                      max.steps=max.steps,
                      weights=weights)

      beta.fit=cbind(fit_L$beta,fit_R$beta)
      coef0.fit=c(fit_L$coef0,fit_R$coef0)
      dev=c(fit_L$deviance,fit_R$deviance)
      lambda=c(fit_L$lambda,fit_R$lambda)

      beta0M=fit_R$beta
      coef0M=fit_R$coef0
      s.list=c(sL,sR)
      aic=c(fit_L$AIC,fit_R$AIC)
      bic=c(fit_L$BIC,fit_R$BIC)
      ebic=c(fit_L$EBIC,fit_R$EBIC)
      while(k<=K.max)
      {
        sM <- round(sL + (sR-sL)*0.618)
        s.list=c(s.list,sM)
        fit_M=bess.glm(x=x,y=y,
                       beta0=beta0M,
                       intercept=coef0M,
                       s=sM,
                       glm.max=glm.max,
                       max.steps=max.steps,
                       weights=weights)
        cat(k,"-th iteration s.left:",sL," s.split:",sM," s.right:",sR,"\\n",sep="")

        beta0M=fit_M$beta
        beta.fit=cbind(beta.fit,beta0M)

        coef0M=fit_M$coef0
        coef0.fit=c(coef0.fit,coef0M)

        dev=c(dev,fit_M$deviance)
        lambda=c(lambda,fit_M$lambda)
        aic=c(aic, fit_M$AIC)
        bic=c(bic, fit_M$BIC)
        ebic=c(ebic, fit_M$EBIC)

        if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev*(sM-sL)) > epsilon &
           abs(fit_R$deviance-fit_M$deviance)/abs(nulldev*(sM-sR)) < epsilon)
        {
          sR <- sM
          fit_R=fit_M
        }else if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev*(sM-sL)) > epsilon &
                 abs(fit_R$deviance-fit_M$deviance)/abs(nulldev*(sM-sR)) > epsilon)
        {
          sL <- sM
          fit_L=fit_M
        }else
        {
          sR=sM
          fit_R=fit_M
          sL=s.min
          fit_L=fit_L1
        }

        if(sR-sL==1) break
        fit_ML=bess.glm(x=x,y=y,
                        beta0=beta0M,
                        intercept=coef0M,
                        s=sM,
                        glm.max=glm.max,
                        max.steps=max.steps,
                        weights=weights)

        fit_MR=bess.glm(x=x,y=y,
                        beta0=beta0M,
                        intercept=coef0M,
                        s=sM,
                        glm.max=glm.max,
                        max.steps=max.steps,
                        weights=weights)
        #if(abs(fit_ML$deviance-fit_M$deviance)/abs(fit_M$deviance) > epsilon &
        # abs(fit_MR$deviance-fit_M$deviance)/abs(fit_M$deviance) < epsilon)
        # {break}
        if(abs(fit_ML$deviance-fit_M$deviance)/abs(nulldev) > epsilon &
           abs(fit_MR$deviance-fit_M$deviance)/abs(nulldev) < epsilon)
        {break}

        k=k+1
      }
    }

    if(method=="sequential")
    {
      nulldev=-2*sum(weights*(y*log(0.5) + (1-y)*log(0.5)))

     # if(abs(dev_L/nulldev)>0.5) dev_L=0
      dev=vector()
      lambda=vector()
      aic=vector()
      bic=vector()
      ebic=vector()
      for(k in 1:length(s.list))
      {
        #cat("select",s.list[k],"variables","\\n")
        if(is.null(factor)){
          if(k == 1){
            fit=bess.glm(x=x,y=y,beta0=beta0,intercept=intercept,
                         s=s.list[k],max.steps=max.steps,glm.max=glm.max,
                         weights=weights)
            beta.fit = matrix(fit$beta)
            coef0.fit = fit$coef0
          }else{
            fit=bess.glm(x=x,y=y,beta0=beta.fit[,k-1,drop=TRUE],intercept=coef0.fit[k-1],
                         s=s.list[k],max.steps=max.steps,glm.max=glm.max,
                         weights=weights)
            beta.fit = cbind(beta.fit,fit$beta)
            coef0.fit = c(coef0.fit,fit$coef0)
          }
        }else{
          if(k == 1){
            fit=gbess.glm(x=x,y=y,Gi=Gi,beta0=beta0,intercept=intercept,
                         s=s.list[k],max.steps=max.steps,glm.max=glm.max,
                         weights=weights)
            beta.fit = matrix(fit$beta)
            coef0.fit = fit$coef0
            s.list[k] = fit$gr_size
          }else{
            fit=gbess.glm(x,y,Gi=Gi,beta0=beta.fit[,k-1,drop=TRUE],intercept=coef0.fit[k-1],
                          s=s.list[k],max.steps=max.steps,glm.max=glm.max,
                          weights=weights)
            beta.fit = cbind(beta.fit,fit$beta)
            coef0.fit = c(coef0.fit,fit$coef0)
            s.list[k] = fit$gr_size
          }
        }
        dev[k]=fit$deviance
        lambda[k]=fit$lambda
        aic[k]=fit$AIC
        bic[k]=fit$BIC
        ebic[k]=fit$EBIC
      }
    }

    beta.fit=sqrt(n)*(beta.fit)/normx
    colnames(beta.fit) = s.list
    rownames(beta.fit) = vn
    coef0.fit = coef0.fit-drop(t(beta.fit)%*%meanx)
    names(coef0.fit) = s.list

    xbest=xs[,which(beta.fit[,ncol(beta.fit),drop=TRUE]!=0)]
    bestmodel=glm(y~xbest, family=binomial, weights=weights)

    if(!setequal(y_names,c(0,1)))
    {
      out=list(family="bess_binomial",method=method,beta=beta.fit,coef0=coef0.fit,s.list=s.list,
               meanx=meanx,normx=normx,nsample=n,bestmodel=bestmodel,
               deviance=dev,nulldeviance=nulldev,AIC=aic,BIC=bic,EBIC=ebic,
               lambda=lambda,y_names=y_names,max.steps=max.steps,factor=factor)
      class(out)="bess"
      return(out)
    }else
    {
      out=list(family="bess_binomial",method=method,beta=beta.fit,coef0=coef0.fit,s.list=s.list,
               meanx=meanx,normx=normx,nsample=n,bestmodel=bestmodel,
               deviance=dev,nulldeviance=nulldev,AIC=aic,BIC=bic,EBIC=ebic,
               lambda=lambda,max.steps=max.steps,factor=factor)
      class(out)="bess"
      return(out)
    }
  }

  if(family=="cox")
  {
    #normalize
    mark=order(y[,1],decreasing = FALSE)
    y=y[mark,]
    weights=weights[mark]
    x=x[mark,]
    beta0=rep(0,p)
    gc()

    if(method=="gsection")
    {
      k = 1
      sL=s.min
      sR=s.max
      beta0R=beta0
      beta0L=beta0

      fit_L1=bess.cox(x,y,
                      beta0=beta0L,
                      s=sL,
                      cox.max=cox.max,
                      max.steps=max.steps,
                      weights=weights)

      nulldev=fit_L1$nulldeviance
      fit_L=fit_L1
      fit_R=bess.cox(x,y,
                     beta0=beta0R,
                     s=sR,
                     cox.max=cox.max,
                     max.steps=max.steps,
                     weights=weights)

      beta.fit=cbind(fit_L$beta,fit_R$beta)
      dev=c(fit_L$deviance,fit_R$deviance)
      lambda=c(fit_L$lambda,fit_R$lambda)

      beta0M=fit_R$beta
      s.list=c(sL,sR)
      aic=c(fit_L$AIC,fit_R$AIC)
      bic=c(fit_L$BIC,fit_R$BIC)
      ebic=c(fit_L$EBIC,fit_R$EBIC)
      while(k<=K.max)
      {
        sM <- round(sL + (sR-sL)*0.618)
        s.list=c(s.list,sM)
        fit_M=bess.cox(x,y,
                       beta0=beta0M,
                       s=sM,
                       cox.max=cox.max,
                       max.steps=max.steps,
                       weights=weights)
        cat(k,"-th iteration s.left:",sL," s.split:",sM," s.right:",sR,"\\n",sep="")

        beta0M=fit_M$beta
        beta.fit=cbind(beta.fit,beta0M)

        dev=c(dev,fit_M$deviance)
        lambda=c(lambda,fit_M$lambda)
        aic=c(aic, fit_M$AIC)
        bic=c(bic, fit_M$BIC)
        ebic=c(ebic, fit_M$EBIC)
        if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev*(sM-sL)) > epsilon &
           abs(fit_R$deviance-fit_M$deviance)/abs(nulldev*(sM-sR)) < epsilon)
        {
          sR <- sM
          fit_R=fit_M
        }else if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev) > epsilon &
                 abs(fit_R$deviance-fit_M$deviance)/abs(nulldev) > epsilon)
        {
          sL <- sM
          fit_L=fit_M
        }else
        {
          sR=sM
          fit_R=fit_M
          sL=s.min
          fit_L=fit_L1
        }

        if(sR-sL==1) break
        fit_ML=bess.cox(x,y,
                        beta0=beta0M,
                        s=sM,
                        cox.max=cox.max,
                        max.steps=max.steps,
                        weights=weights)

        fit_MR=bess.cox(x,y,
                        beta0=beta0M,
                        s=sM,
                        cox.max=cox.max,
                        max.steps=max.steps,
                        weights=weights)
        if(abs(fit_ML$deviance-fit_M$deviance)/abs(fit_M$deviance) > epsilon &
           abs(fit_MR$deviance-fit_M$deviance)/abs(fit_M$deviance) < epsilon)
        {break}

        k=k+1
      }
    }

    if(method=="sequential")
    {
      dev=vector()
      lambda=vector()
      aic=vector()
      bic=vector()
      ebic=vector()
      for(k in 1:length(s.list))
      {
        #cat("select",s.list[k],"variables","\\n")
        if(is.null(factor)){
          if(k == 1){
            fit=bess.cox(x=x,y=y,beta0=beta0,
                         s=s.list[k],cox.max=cox.max,
                         max.steps=max.steps,weights=weights)
            beta.fit = matrix(fit$beta)
          }else{
            fit=bess.cox(x=x,y=y,beta0=beta.fit[,k-1,drop=TRUE],
                         s=s.list[k],cox.max=cox.max,
                         max.steps=max.steps,weights=weights)
            beta.fit = cbind(beta.fit,fit$beta)
          }
        }else{
          if(k == 1){
            fit=gbess.cox(x=x,y=y,Gi=Gi,beta0=beta0,
                          s=s.list[k],cox.max=cox.max,
                          max.steps=max.steps,weights=weights)
            beta.fit = matrix(fit$beta)
            s.list[k] = fit$gr_size
          }else{
            fit=gbess.cox(x,y,Gi=Gi,beta0=beta.fit[,k-1,drop=TRUE],
                          s=s.list[k],cox.max=cox.max,
                          max.steps=max.steps,weights=weights)
            beta.fit = cbind(beta.fit,fit$beta)
            s.list[k] = fit$gr_size
          }
        }
        dev[k]=fit$deviance
        lambda[k]=fit$lambda
        aic[k]=fit$AIC
        bic[k]=fit$BIC
        ebic[k]=fit$EBIC
      }
    }
    beta.fit=sqrt(n)*(beta.fit)/normx
    colnames(beta.fit) = s.list
    rownames(beta.fit) = vn

    xbest=xs[,which(beta.fit[,ncol(beta.fit),drop=TRUE]!=0)]
    bestmodel=coxph(Surv(y[,1],y[,2])~xbest, iter.max=cox.max, weights=weights)
    nulldev = bestmodel$loglik[1]

    out=list(family="bess_cox",method=method,beta=beta.fit,s.list=s.list,meanx=meanx,
             normx=normx,nsample=n,bestmodel=bestmodel,
             deviance=dev,nulldeviance=nulldev,AIC=aic,BIC=bic,EBIC=ebic,
             lambda=lambda,max.steps=max.steps,factor=factor)
    class(out)="bess"
    return(out)
  }

}
