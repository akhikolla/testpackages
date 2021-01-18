  #' linear model based on CLV 
  #' 
  #' prediction of a response variable, y, based on clusters of predictors variables, X.
  #' boosted-liked procedure for identifying groups of predictors, 
  #' and their associated latent component,
  #' well correlated with the actual residuals of response variable, y.
  #' sparsity is allowed using the strategy options ("sparselv" or "kplusone") and the rho parameter. 
  
  #' @param X : The matrix of the predictors, the variables to be clustered
  #' @param y : The response variable (usually numeric)
  #'            If y is binary factor, indicator variable (0/1) is generated. A Bayes rule is used to compute class probabilities.\cr
  #'            Performance criteria is RMSE for numerical variable; RMSE and error rate for binary factor. \cr
  #' @param method : The criterion to be use in the cluster analysis.  \cr
  #'        1 or "directional" : the squared covariance is used as a measure of proximity (directional groups). \cr    
  #'        2 or "local"       : the covariance is used as a measure of proximity (local groups)
  #' @param sX : TRUE/FALSE, i.e. standardization or not of the columns X (TRUE by default)
  #' @param shrinkp : shrinkage  paramater used in the boosting (max : 1, 0.5 by default). \cr
  #'         If shrinkp is a vector of positive values greater than 0, and lower or equal to 1, the outputs are given for each value.
  #' @param strategy : "none" (by default), or "kplusone" (an additional cluster for the unclassifiable variables),
  #'        or "sparselv" (zero loadings for the unclassifiable variables)
  #' @param rho : a threshold of correlation between 0 and 1 (used in "kplusone" or "sparselv" strategy, 0.3 by default)
  #' @param validation TRUE/FALSE i.e. using a test set or not. By default no validation
  #' @param id.test : if validation==TRUE, the number of the observations used as test set
  #' @param maxiter : the maximum number of components extracted (100 by default)
  #' @param threshold : used in a stopping rule, when the relative calibration errors sum of squares stabilizes (10e-6 by default)
  #' 
  #' @return \item{Group}{a list of the groups of variables X in order of the first time extracted.}
  #' @return \item{Comp}{ a list of the latent components associated with the groups of X variables extracted.}
  #' @return \item{Load}{ a list for the loadings of the X variables in the latent component.}
  #' @return \item{Alpha}{ a list of the regression coefficients to be applied to the latent components. \cr
  #'         The coefficients are aggregated when the same latent component is extracted several times during the iterative steps.}
  #' @return \item{Beta}{ a list of the beta coefficients to be applied to the pretreated predictors. \cr
  #'         For a model with the A first latent components, the A first elements of the list must be added together.}  
  #' @return \item{GroupImp}{ Group Importance i.e. the decrease of the residuals' variance provided by the CLV components in the model.}
  #' @return \item{RMSE.cal}{ the root mean square error for the calibration set, at each step of the procedure.}
  #' @return \item{ERRrate.cal and rocAUC.cal}{ when y is a binary factor, the classification rate and the AUC for ROC, on the bassis of the calibration set, at each step of the procedure.}
  #' @return \item{RMSE.val}{ as RMSE.cal but for the test set, if provided.}
  #' @return \item{ERRrate.val and rocAUC.val}{as for calibration set but for the test set, if provided.}
  #' 
  #' @seealso CLV, CLV_kmeans
  #'
  #' @export
  #' 
lm_CLV = function(X,y,method="directional",sX=TRUE,shrinkp=0.5, strategy="none",rho=0.3, validation=FALSE,id.test=NULL,maxiter=100, threshold=10e-6 )
{  
  
  if(is.null(method)) stop('parameter method should be =1/"directional" or =2/"local"')
  if (method=="directional") method=1
  if (method=="local") method=2
  
  n = dim(X)[1]
  p = dim(X)[2]
  gpmax=p
  #maxiter=maxiter/shrinkp
  
  yisfact=FALSE
  yfact=NULL
  ERRrate.cal=NULL
  ERRrate.val=NULL
  rocAUC.cal=NULL
  rocAUC.val=NULL
  if (is.factor(y)) {
    yfact=y
    y=as.numeric(yfact)-1
    yisfact=TRUE
    n1=sum(y==1)
    n0=sum(y==0)
  }
  
  if (!validation) {
    Xscale=scale(X,scale=sX)
    scale.mean<-apply(X,2,mean)
    scale.scale<-apply(X,2,sd)
    meany=mean(y)
    ycenter=scale(y,scale=FALSE)
  }
  if (validation) {
    test<-id.test
    train<-setdiff(1:n,test)
    meany=mean(y[train])
    scale.mean<-apply(X[train,],2,mean)
    Xtest=X[test,]
    ytest=y[test]
    if (length(test)==1) Xtest=t(Xtest)
    if (sX==TRUE) {
      scale.scale<-apply(X[train,],2,sd)
    } else {
      scale.scale=rep(1,ncol(X))
    }
    X=X[train,]
    y=y[train]
    n<-length(train)
    Xscale=scale(X,scale=sX)
    ycenter=scale(y,scale=FALSE)
    ERR.val<-matrix(ytest-rep(meany,length(ytest)),length(id.test),1)
  }
  ERR.cal<-matrix(y-rep(mean(y),length(y)),n,1)
 
  # CLV  -----------------------------------------------
  # performed on X, one time only
  # ascendant hierarchy, with consolidation
  #start.time <- Sys.time()
  res.clv_X=CLV(X,sX=sX,method=method,nmax=gpmax)  # avec consolidation 
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time;          print(time.taken)
  
# ==============================================================  
# for each value of the shrinkage parameter 
  resul=list()
  for (valnu in 1:length(shrinkp))  {
    Group = list()
    Load=list()
    Comp =list()
    Beta=list() #matrix(0,maxiter,ncol(X))
    Alpha=list()
 
 # ------------------------------------------------------------
 # iterative search for groups of X variables (and associated LV) well correlated with y
                        # yt=rep(mean(y),length(y))
  et=ycenter
  iter=0
                         # stockyt=yt
  stocket=et
  Xscale=scale(X,scale=sX)
 
  while (iter<maxiter) {
      iter=iter+1
                               #print(iter)
      tabgpy=matrix(0,p,gpmax)
      idgpy=rep(0,gpmax)
      compgpy=matrix(0,n,gpmax)
      relate = list()
      for (k in 1:gpmax) {
        #gp=res.clv_X[[k]]$clusters[1,]
        gp=get_partition(res.clv_X,k)
        comp=get_comp(res.clv_X,k)
        if (method==1) gpy=which.max(abs(cor(comp,et)))
        if (method==2) gpy=which.max(cor(comp,et))
        idgpy[k]=gpy
        compgpy[,k]=comp[,gpy]
        tabgpy[which(gp==gpy),k]=1
        relate[[k]]=which(gp==gpy)
      }
  
      # ------------------------------------------------------------------------
      # identify the levels for which the X variables associated with y change
      # ------------------------------------------------------------------------
      nivparti=NULL
      for (k in 1:(length(relate)-1)) {
        respons=length(setdiff(relate[[k]],relate[[k+1]]))+length(setdiff(relate[[k+1]],relate[[k]]))
        if (respons>0) nivparti=c(nivparti,k) 
      }
      memoiter=cbind(nivparti,idgpy[nivparti])
      memorelate=relate[nivparti]
      memocomp=compgpy[,nivparti]
                                            
  

      # ------------------------------------------------------------------------
      # updating of groups and latent components if strategy <>"none
      # only the group "gpy" at each level of "nivparti", and only one iteration (for simplicity)
      # ------------------------------------------------------------------------
      if (strategy=="sparselv"){
        for (kk in 1:length(nivparti)) {
             res=consol_calcul_s(method=method,X=Xscale,EXTr=0,Xr=NULL,EXTu=0,Xu=NULL,
                              ind=memorelate[[kk]],max.iter=20, eps = 0.001,rlevel=rho)
             memocomp[,kk]<-as.matrix(res$comp)
             sloading<-as.matrix(res$loading)
             memorelate[[kk]]=memorelate[[kk]][which(sloading!=0)]
         }
      }
      if (strategy=="kplusone"){
        for (kk in 1:length(nivparti)) {
          Xkk=as.matrix(Xscale[,memorelate[[kk]]])                          
          gptemp=consol_affect_k(method=method,X=Xkk,Xr=NULL,Xu=NULL,EXTr=0,EXTu=0,comp=memocomp[,kk], a=NULL, u=NULL,rlevel=rho)
          memorelate[[kk]]=memorelate[[kk]][which(gptemp==1)]
          res<-consol_calcul(method=method,X=Xscale,EXTr=0,Xr=NULL,EXTu=0,Xu=NULL,ind=memorelate[[kk]]) 
          memocomp[,kk]=res$comp
        }
      }
      # ------------------------------------------------------------------------
      # choice of one cluster according to the modified Kaiser Guttman criterion
      # ------------------------------------------------------------------------
      rescrit=NULL
      for(j in 1:length(nivparti)){
        indiy=memorelate[[j]]
                                           #print(indiy)
        py=length(indiy)
        if (py==0) {break}
        # let the X variables associated with y : is this group unidimensionnal ?
        Xj=as.matrix(X[,indiy])
        Xjsc=scale(Xj)
        svdXjsc=svd(Xjsc/sqrt(n-1))
        # KG
        l1=svdXjsc$d[1]^2
        l2=svdXjsc$d[2]^2
        seuill=1+ ( 2* sqrt((py-1)/(n-1)))
        KGX=(l1>seuill)&(l2<seuill)
        rescrit=rbind(rescrit,c(aclevel=nivparti[j],py=py,l1=l1,l2=l2,seuill=seuill,KGX=KGX))
      }
      #tabgpyniv=reschoice$tabgpyniv
      rescrit=as.data.frame(rescrit) 
      change=which(diff(rescrit$KGX)==1)
      change=change[length(change)]+1     # with this option, selection of the first step the groups remain unidimensional
      nf=rescrit[change,1];
      if (length(nf)==0) {nf=nivparti[length(nivparti)]}   # if no solution nf found, choice of the last solution

      # ------------------------------------------------------------------------
      # extraction of the Group and Comp, Evaluation of Load at this iteration/step
      # ------------------------------------------------------------------------
      #if (!is.na(nf)) {
          ind=memorelate[[which(nivparti==nf)]]
          indname=colnames(X)[ind]
          group=indname
          Xiter=as.matrix(Xscale[,ind])
          if (method==1) {
            #peig=powerEigen(crossprod(Xiter))
            #load=as.matrix(peig$vectors)
            reseig=eigen(crossprod(Xiter))
            load=as.matrix(reseig$vectors[,1])
            sig=sign(load[which.max(load)])
            load=load*sig
          }
          if (method==2) {
            load=matrix(1/length(ind),length(ind),1)
          }
          rownames(load)=colnames(Xiter)
          comp=Xiter%*%load      # idem...but the sign as : as.matrix(memocomp[,which(nivparti==nf)]) 
          #res=extraction(nivparti,tabgpyniv,nf,res.clvX,strategy,Xscale)
          Group[[iter]]=group
          Load[[iter]]=load
          Comp[[iter]]=comp
      # } else {
      #     break
      # }
      
      
      # ---------------------------------------------------------------
      # boosted-like step : beta, prediction and deflation
      # --------------------------------------------------------------
          alpha=shrinkp[valnu]*(cov(et,comp)/var(comp)) 
          Alpha[[iter]]=c(alpha)
          # ouind=apply(sapply(Group,FUN="==",colnames(X)),1,sum)
          # ind=which(ouind==1)
          Beta[[iter]]=rep(0,ncol(X))
          Beta[[iter]][ind]=as.numeric(alpha)*load
          yiter=comp%*%alpha
          et= as.vector( et - yiter )   
          stocket=cbind(stocket,et);
                                # yt=yt+yiter
                                # stockyt=cbind(stockyt,yt)
                                # print(stockyt)
                                #memo=memoiter[1:which(nivparti==nf),]
       
        # stopping rule                                         
        # if relative variance of et < threshold (10e-6 by default)
        critere=(var(stocket[,iter])-var(stocket[,iter+1]))  /var(y)
                                         #print(critere)
         if (critere<threshold) {
              if (critere<0) print("Warnings regarding !")
              break
         }
       
    
  } # end of iter
  nbiter=iter
                         print(paste("nb iterations: ",nbiter))
  
  ##############################################################
  # if the same group has been extracted several times
  # the outputs are simplified, 
  # each group is listed only one time (by order of the first apparition)
  # the coefficients Alpha and Beta, are aggregated 
  # a Group Importance criterion is evaluated
   uGroup=unique(Group)
   uLoad=list()
   uComp=list()
   uAlpha=list()
   uBeta=list()
   w=-diff(apply(stocket,2,var)) 
   #w=w/c(var(y))
   
   resdetail=matrix(0,nbiter,length(uGroup))                       
   for (g in 1:length(uGroup)) {
   cas=0
     for (iter in 1:length(Group)) {
          if(identical(uGroup[[g]],Group[[iter]])) {
            cas=cas+1
            resdetail[iter,g]=w[iter]
            if (cas==1) {
              uLoad[[g]]=Load[[iter]]
              uComp[[g]]=Comp[[iter]]
              uAlpha[[g]]=Alpha[[iter]]
              uBeta[[g]]=Beta[[iter]]
            } else {
              uAlpha[[g]]=uAlpha[[g]]+Alpha[[iter]]
              uBeta[[g]]=uBeta[[g]]+Beta[[iter]]
            }
         }
     }
   }
    GroupImp=apply(resdetail,2,sum)
    # re_order the group according to their importance
    ordr=order(GroupImp,decreasing=TRUE)
    GroupImp=GroupImp[ordr]
    uGroup <- uGroup[ordr]
    uLoad <- uLoad[ordr]
    uComp <- uComp[ordr]
    uAlpha <- uAlpha[ordr]
    uBeta <- uBeta[ordr]
    
  # ##########################################
  # final section : prediction
  RMSE.cal=NULL
  if (yisfact) { ERRrate.cal=NULL; rocAUC.cal=NULL} 
  # RMSP.cal=NULL  
  if (validation)    {
    RMSE.val=NULL
    if (yisfact) { ERRrate.val=NULL; rocAUC.val=NULL} 
  # RMSP.val=NULL
    Xstest<-scale(Xtest,center=scale.mean,scale=scale.scale)
  }
  # step "cst model"
  ypred=rep(mean(y),length(y))
  RMSE.cal=sqrt(mean((y-ypred)^2))
  if (yisfact) {
    m1=mean(ypred[y==1]);  s1=sd(ypred[y==1]);
    m0=mean(ypred[y==0]);  s0=sd(ypred[y==0]);
    tab=table(y,Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$yclass)
    ERRrate.cal=(length(y)-sum(diag(tab)))/length(y)
    W=wilcox.test(Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$pclass~ y,exact=FALSE)$statistic
    rocAUC.cal=1-(W/(sum(y==1)*sum(y==0)))
  }
  # RMSP.cal=c(sqrt(mean(scale(ypred,scale=FALSE)^2)))
  if (validation) {
    ypredtest=rep(meany,length(ytest))
    RMSE.val=c(sqrt(mean((ytest-ypredtest)^2)))
    if (yisfact) {
      tab=table(ytest,Bayes_classif(ypredtest,n1,m1,s1,n0,m0,s0)$yclass)
      ERRrate.val=(length(ytest)-sum(diag(tab)))/length(ytest)
      W=wilcox.test(Bayes_classif(ypredtest,n1,m1,s1,n0,m0,s0)$pclass~ ytest,exact=FALSE)$statistic
      rocAUC.val=1-(W/(sum(ytest==1)*sum(ytest==0)))
    }
  # RMSP.val=c(sqrt(mean(scale(ypredtest,scale=FALSE)^2)))
  }
  # step by step
  for (iter in 1 :nbiter ) {
     ypred=ypred+Xscale%*%Beta[[iter]]
     #ERR.cal=cbind(ERR.cal,y-ypred)
     RMSE.cal=c(RMSE.cal,sqrt(mean((y-ypred)^2)))
     if (yisfact) {
       m1=mean(ypred[y==1]);  s1=sd(ypred[y==1]);
       m0=mean(ypred[y==0]);  s0=sd(ypred[y==0]);
       tab=table(y,Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$yclass)
       ERRrate.cal=c(ERRrate.cal,(length(y)-sum(diag(tab)))/length(y))
       W=wilcox.test(Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$pclass~ y,exact=FALSE)$statistic
       rocAUC.cal=c(rocAUC.cal,1-(W/(sum(y==1)*sum(y==0))))
     }
     #RMSP.cal=c(RMSP.cal,sqrt(mean(scale(ypred,scale=FALSE)^2)))
     if (validation){
          ypredtest=ypredtest+Xstest%*%Beta[[iter]]
          #ERR.val=cbind(ERR.val,ytest-ypredtest)
          RMSE.val=c(RMSE.val,sqrt(mean((ytest-ypredtest)^2)))
          if (yisfact) {
            tab=table(ytest,Bayes_classif(ypredtest,n1,m1,s1,n0,m0,s0)$yclass)
            ERRrate.val=c(ERRrate.val,(length(ytest)-sum(diag(tab)))/length(ytest))
            W=wilcox.test(Bayes_classif(ypredtest,n1,m1,s1,n0,m0,s0)$pclass~ ytest,exact=FALSE)$statistic
            rocAUC.val=c(rocAUC.val,1-(W/(sum(ytest==1)*sum(ytest==0))))
          }
          #RMSP.val=c(RMSP.val,sqrt(mean(scale(ypredtest,scale=FALSE)^2)))
     } 
   }
                  
    if (!validation)  resul[[valnu]]=list(shrinkp=shrinkp[valnu],Group=uGroup,Load=uLoad,Comp=uComp,Alpha=uAlpha,Cst=meany,Beta=uBeta,GroupImp=GroupImp,RMSE.cal=RMSE.cal,ERRrate.cal=ERRrate.cal,rocAUC.cal=rocAUC.cal,resdetail=resdetail,sX=sX,X=X,yfact=yfact)
    if (validation)   resul[[valnu]]=list(shrinkp=shrinkp[valnu],Group=uGroup,Load=uLoad,Comp=uComp,Alpha=uAlpha,Cst=meany,Beta=uBeta,GroupImp=GroupImp,RMSE.cal=RMSE.cal,ERRrate.cal=ERRrate.cal,rocAUC.cal=rocAUC.cal,RMSE.val=RMSE.val,ERRrate.val=ERRrate.val,rocAUC.val=rocAUC.val, resdetail=resdetail,sX=sX,X=X,yfact=yfact)
 
  }
  # end of for the valnu values
  #=======================================================================
  
  class(resul) = "lmclv"
  return(resul)
}


