

#####################################################
#####  2-Parameter Logistic Model with L2-Norm  #####
#####################################################

P2L2 = function(x, y, dose0, tp=0.3, mlambda=1, nlambda=50, rlambda=NULL, wldose=NULL, nfolds=length(y), foldid=NULL, keep.beta=FALSE, thresh=1e-6, maxit=1e+3, threshP=1e-5, threshB=100) {

  storage.mode(y)="double"; N0=length(y)
  x1=cbind(1.0, log(x)); p1=ncol(x1)
  ltp=log(tp/(1.0-tp))
  ldose0=log(dose0)
  nd=length(dose0)


  if (all(y==0) | all(y==1)) {

    if (all(y==0)) {
      betai=c(log(0.001/0.999),0.001)
    } else {
      betai=c(log(0.999/0.001),0.001)
    }

    probi=1/(1+exp(betai[1]+betai[2]*ldose0)^{-1})
    temi=abs(probi-tp); id=which.min(temi)

    return(list(Beta=matrix(betai,ncol=1), prob=probi, dose.closest=id))

  } else {

    ### Weights - MLE ###
    if (is.null(wldose)) {
      mlej=P2L2(x,y,dose0,tp,mlambda=0,nlambda=1,rlambda,wldose=rep(1, nd),nfolds=1,foldid,keep.beta=FALSE,thresh,maxit,threshP,threshB)
      betaj=mlej$Beta

      probj0=1/(1+exp(betaj[1,1]+betaj[2,1]*ldose0)^{-1})
      wldose=pmax(0.001,(betaj[1,1]+betaj[2,1]*ldose0-ltp)^2)

      wldose=1.0/wldose
      wldose=wldose/sum(wldose)*nd
    }


    #####  Run  #####
    x1i=x1; p1i=p1
    indexi=which(apply(x1i,2,sd)==0)[-1]
    if (length(indexi)>0) {
      x1i=as.matrix(x1i[,-indexi],nrow=N0); p1i=p1i-length(indexi)
    }

    ### Lambda path ###
    if (is.null(rlambda)) {
      rlambda=ifelse(N0>p1i, 0.0001, 0.01)
    }
    if (nlambda==1 | mlambda==0) {
      lambda=mlambda; nlambda=1
    } else {
      lambda=mlambda*(rlambda)^(c(0:(nlambda-1))/(nlambda-1))
    }


    out=P2L2C(x1i, y, ldose0, ltp, lambda, wldose, thresh, maxit, threshP, threshB)

    nlambdai=out$nlambda ## number of lambdas
    if (nlambdai==0)
      return(NULL)
    lambdai=out$lambda[1:nlambdai]

    if (length(indexi)>0) {
      betai=matrix(0,nrow=p1,ncol=nlambdai)
      betai[-indexi,]=out$Beta
      out$Beta=betai

      betai=matrix(0,nrow=p1,ncol=nlambdai)
      betai[-indexi,]=out$BetaSTD
      out$BetaSTD=betai
    }

    out$Beta=matrix(out$Beta[, 1:nlambdai], ncol=nlambdai)
    out$BetaSTD=matrix(out$BetaSTD[, 1:nlambdai], ncol=nlambdai)
    # out$nzero=apply(out$Beta!=0, 2, sum)
    out$flag=out$flag[1:nlambdai]
    out$LL=out$LL[1:nlambdai]

    if ((nfolds==1 & is.null(foldid)) | nlambda==1) {
      fit=data.frame(lambda=lambdai, pDev=(out$ll0-out$LL)/out$ll0)

      if (nlambda==1) {
        betai=c(out$Beta)
        probi=1/(1+exp(betai[1]+betai[2]*ldose0)^{-1})
        temi=abs(probi-tp); id=which.min(temi)

        return(list(Beta=out$Beta, fit=fit, flag=out$flag, prob=probi, dose.closest=id))
      }

      return(list(Beta=out$Beta, fit=fit, flag=out$flag))
      # return(list(Beta=out$BetaSTD, fit=fit, flag=out$flag))

    } else {

      ###  Split data for cross-validation  ###
      if (is.null(foldid)) {
        foldid=sample(rep(seq(nfolds), length=N0))
      } else {
        nfolds=max(foldid)
      }
      tb=table(foldid)
      N0i=numeric(nfolds); Nf=numeric(nfolds)
      for (i in 1:nfolds) {
        N0i[i]=sum(tb[-i]); Nf[i]=tb[i]
      }
      weighti=as.vector(tapply(rep(1,N0), foldid, sum))


      #####  Cross-validation estimates  #####
      # ypred=matrix(0,nrow=N0,ncol=nlambdai)

      outi=list(); cvRSS=matrix(NA, nrow=nfolds, ncol=nlambdai); i=1
      for (i in 1:nfolds) {
        temid=(foldid==i)

        if (any(y[!temid]==0) & any(y[!temid]==1)) {
          x1i=matrix(x1[!temid,,drop=FALSE],nrow=N0i[i]); x1j=matrix(x1[temid,,drop=FALSE],nrow=Nf[i]); p1i=p1

          indexi=which(apply(x1i,2,sd)==0)[-1]
          if (length(indexi)>0) {
            x1i=matrix(x1i[,-indexi],nrow=N0i[i]); x1j=matrix(x1j[,-indexi],nrow=Nf[i])
            p1i=p1i-length(indexi)
          }

          outi[[i]]=cvP2L2C(x1i, y[!temid], x1j, y[temid], ldose0, ltp, lambdai, wldose, thresh, maxit, threshP, threshB)

          cvRSS[i, 1:outi[[i]]$nlambda]=2*(0-outi[[i]]$LLF)*Nf[i] ## for ith fold
        }
      }

      cvRSS=matrix(cvRSS[, 1:nlambdai], ncol=nlambdai)
      cvraw=cvRSS/weighti; nfoldi=apply(!is.na(cvraw), 2, sum); #rm(cvRSS) #
      cvm=apply(cvraw, 2, weighted.mean, w=weighti, na.rm=TRUE)
      cvse=sqrt(apply(sweep(cvraw, 2, cvm, "-")^2, 2, weighted.mean, w=weighti, na.rm=TRUE)/(nfoldi-1))

      indexi=which.min(cvm)
      indexij=which(cvm<=(cvm[indexi]+cvse[indexi]))[1]
      temi=rep("", nlambdai)
      temi[indexi]="*";#temi[indexij]=ifelse(temi[indexij]=="", "*", "***")
      #temCV=data.frame(lambda=lambdai, cvm=cvm, cvse=cvse, nzero=out$nzero, index=temi,stringsAsFactors=FALSE)
      temCV=data.frame(lambda=lambdai, pDev=(out$ll0-out$LL)/out$ll0, cvm=cvm, cvse=cvse, index=temi, stringsAsFactors=FALSE)


      betai=c(out$Beta[,indexi])
      probi=1/(1+exp(betai[1]+betai[2]*ldose0)^{-1})
      temi=abs(probi-tp); id=which.min(temi)

      rm(outi)
      if (!keep.beta) {
        # lambda.1se=lambdai[indexij]
        return(list(Beta=matrix(out$Beta[, indexi],ncol=1), fit=temCV, lambda.min=lambdai[indexi], flag=out$flag, prob=probi, dose.closest=id))
      } else {
        return(list(Beta=out$Beta, fit=temCV, lambda.min=lambdai[indexi], flag=out$flag, prob=probi, dose.closest=id))
      }

    } # folder

  } # rCRM

}


