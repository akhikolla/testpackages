#' @title Computationally efficient conservative estimate
#'
#' @description Computes conservative estimates with two step GANMC procedure for a Gaussian vector. The probability is approximated with a biased low dimensional estimator and the bias is corrected with a MC estimator.
#'
#' @param alpha probability of conservative estimate.
#' @param pred list containing mean vector (pred$mean) and covariance matrix (pred$cov).
#' @param design a matrix of size \code{length(pred$mean)}x(input space dimension) that contains the design where \code{pred$mean} was computed.
#' @param threshold threshold, real number.
#' @param pn coverage probability function, vector of the same length as pred$mean (if not specified it is computed).
#' @param type type of excursion: ">" for excursion above threshold or "<" for below.
#' @param verb level of verbosity, integer from 1--7.
#' @param lightReturn boolean for light return.
#' @param algo choice of algorithm for computing probabilities ("GANMC", "GMC").
#' @return A list containing the conservative estimate (\code{set}), the Vorob'ev level (\code{lvs}). If \code{lightReturn=FALSE}, it also returns the actual probability of the set (\code{proba}) and the variance of this estimate (\code{vars}).
#' @examples
#' if (!requireNamespace("DiceKriging", quietly = TRUE)) {
#' stop("DiceKriging needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # Compute conservative estimate of excursion set of testfun below threshold
#' # Initialize
#' testfun<-function(x){return(((3*x^2+7*x-3)*exp(-1*(x)^2)*cos(5*pi*x^2)-1.2*x^2))}
#' mDet<- 1500
#'
#' # Uniform design points
#' set.seed(42)
#' doe<-runif(n = 8)
#' res<-testfun(doe)
#' threshold<-0
#' # create km
#' smallKm <- DiceKriging::km(design = matrix(doe,ncol=1),
#' response = res,covtype = "matern5_2",coef.trend = -1,coef.cov = c(0.05),coef.var = 1.1)
#' # prediction at newdata
#' newdata<-data.frame(matrix(seq(0,1,,mDet),ncol=1)); colnames(newdata)<-colnames(smallKm@X)
#' pred<-DiceKriging::predict.km(object = smallKm,newdata = newdata,type = "UK",cov.compute = TRUE)
#'
#' \dontrun{
#' # Plot (optional)
#' plot(seq(0,1,,mDet),pred$mean,type='l')
#' points(doe,res,pch=10)
#' abline(h = threshold)
#' lines(seq(0,1,,mDet),pred$mean+pred$sd,lty=2,col=1)
#' lines(seq(0,1,,mDet),pred$mean-pred$sd,lty=2,col=1)
#' }
#' # Compute the coverage function
#' pn<-pnorm((threshold-pred$mean)/pred$sd)
#'
#' \dontrun{
#' pred$cov <- pred$cov + 1e-7*diag(nrow = nrow(pred$cov),ncol = ncol(pred$cov))
#' CE<-conservativeEstimate(alpha = 0.95,pred = pred,design = as.matrix(newdata),
#' threshold = threshold,type = "<",verb=1, pn=pn,algo = "ANMC")
#' points(newdata[CE$set,],rep(-0.1,mDet)[CE$set],col=4,pch="-",cex=2)
#' }
#' @references Azzimonti, D. and Ginsbourger, D. (2018). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Journal of Computational and Graphical Statistics, 27(2), 255-267. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#'
#' Bolin, D. and Lindgren, F. (2015). Excursion and contour uncertainty regions for latent Gaussian models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(1):85--106.
#' @export
conservativeEstimate<-function(alpha=0.95,pred,design,threshold,pn=NULL,type=">",verb=1,lightReturn=T,algo="GANMC"){
  if(is.null(pn)){
    sds<-sqrt(diag(pred$cov))
    if(type==">"){
      pn <- pnorm((pred$mean -threshold)/sds)
    }else{
      pn <- pnorm((threshold-pred$mean)/sds)
    }
  }
  sortPn<-sort(pn,index.return=T,decreasing = T)

  productPn<-rep(-1,1e4)
  productPn[1]<-sortPn$x[1]
  indMaxSort<-1
  indxAlphaProd<-1
  while(sortPn$x[indMaxSort]>alpha){
    indMaxSort<-indMaxSort+1
    productPn[indMaxSort]<-prod(sortPn$x[1:indMaxSort])
    if(productPn[indMaxSort]>alpha)
      indxAlphaProd<-indxAlphaProd+1
  }

  if(verb){
    cat("################################### \n")
    cat("Conservative Estimates        \n")
    if(verb>1)
      cat("Alpha: ",alpha, ", type: ",type,", algorithm: ",algo,"\n")
    cat("################################### \n\n")
    cat("Creating big covariance matrix (size:",indMaxSort,")...")
  }
  indx<-sortPn$ix[1:indMaxSort]
  Sigma<-matrix(0,nrow = indMaxSort,ncol = indMaxSort)
  for(k1 in seq(indMaxSort))
    for(k2 in seq(indMaxSort)){
      Sigma[k1,k2]<-pred$cov[indx[k1],indx[k2]]
    }
  if(verb)
    cat("Done.\n Memory currently used: ",sum(gc()[,2]), " MB\n Initializing other variables...")
  E<-matrix(design[indx,],ncol=ncol(design))
  mu<-pred$mean[indx]

  leftIndx<-indxAlphaProd
  rightIndx<-indMaxSort

  meth=4
  sizeMaxGenz<-500

  if(verb)
    cat("Done.\nComputing first step...\n")
  # if indMaxSort is less than 100 there's no need to use MCgrf
  if(indMaxSort<sizeMaxGenz){
    if(type==">"){
      a<-threshold
      b<-Inf
    }else{
      a<--Inf
      b<-threshold
    }
    if(verb>1)
      cat(" Size of excursion (",indMaxSort,") less than",sizeMaxGenz,", use Genz \n")
    ProbaRight<-pmvnorm(lower = a,upper = b,mean = mu[1:rightIndx],sigma = Sigma[1:rightIndx,1:rightIndx])
    ProbaLeft<-pmvnorm(lower = a,upper = b,mean = mu[1:leftIndx],sigma = Sigma[1:leftIndx,1:leftIndx])
    varRight<-attr(ProbaRight,"error")
    varLeft<-attr(ProbaLeft,"error")
  }
  else{ # Otherwise we need to..
    precMC<-2.5e3
    timeGANMC<-10
    qq<-c(min(10,leftIndx),min(20,leftIndx))
    if(verb>1)
      cat(" indMaxSort (",indMaxSort,") greater than",sizeMaxGenz,", use MC/ANMC \n")
    if(type==">"){
      if(verb>1){
        cat("  Type: > \n")
        cat("  Compute probaL... \n")
        cat("#####################\n")
      }
      probaL<-ProbaMin(cBdg=timeGANMC,q = qq,E = E[1:leftIndx,],threshold = threshold,mu = mu[1:leftIndx],Sigma = Sigma[1:leftIndx,1:leftIndx],pn = NULL,lightReturn=T,method=meth,verb=(verb-1),Algo=algo)
      resL<-1-probaL$probability
      if(verb>1){
        cat("#####################\n")
        cat("   Done.\n  Compute probaR... \n")
        cat("#####################\n")
      }
      probaR<-ProbaMin(cBdg=timeGANMC,q = qq,E = E[1:rightIndx,],threshold = threshold,mu = mu[1:rightIndx],Sigma = Sigma[1:rightIndx,1:rightIndx],pn = NULL,lightReturn=T,method=meth,verb=(verb-1),Algo=algo)
      resR<-1-probaR$probability
      if(verb>1){
        cat("#####################\n")
        cat("   Done.\n")
      }
    }else{
      if(verb>1){
        cat("  Type: < \n")
        cat("  Compute probaL...\n")
        cat("#####################\n")
      }
      pn1L<- sortPn$x[1:leftIndx]#pnorm((threshold-mu[1:leftIndx])/sqrt(diag(Sigma[1:leftIndx,1:leftIndx])))
      probaL<-ProbaMax(cBdg=timeGANMC,q = qq,E = E[1:leftIndx,],threshold = threshold,mu = mu[1:leftIndx],Sigma = Sigma[1:leftIndx,1:leftIndx],pn = pn1L,lightReturn=T,method=meth,verb=(verb-1),Algo=algo)
      resL<-1-probaL$probability
      if(verb>1){
        cat("#####################\n")
        cat("Done.\n  Compute probaR... \n")
        cat("#####################\n")
      }
      pn1R<-sortPn$x[1:rightIndx] #pnorm((threshold-mu[1:rightIndx])/sqrt(diag(Sigma[1:rightIndx,1:rightIndx])))
      probaR<-ProbaMax(cBdg=timeGANMC,q = qq,E = E[1:rightIndx,],threshold = threshold,mu = mu[1:rightIndx],Sigma = Sigma[1:rightIndx,1:rightIndx],pn = pn1R,lightReturn=T,method=meth,verb=(verb-1),Algo=algo)
      resR<-1-probaR$probability
      if(verb>1){
        cat("#####################\n")
        cat("Done.\n")
      }
    }
    ProbaLeft<-as.numeric(resL)
    ProbaRight<-as.numeric(resR)
    varLeft<- probaL$variance
    varRight<- probaR$variance
  }
  if(verb)
    cat("Done. (Memory currently used: ",sum(gc()[,2])," MB)\n")
  if(verb){
    cat("\nIteration 0 finished.")
    if(verb>1){
      if(indMaxSort<sizeMaxGenz){
        cat("leftIndx:",leftIndx,",ProbaLeft: ",ProbaLeft,", rightInx: ",rightIndx,",ProbaRight: ",ProbaRight," (Genz)")
      }else{
        cat("leftIndx:",leftIndx,",ProbaLeft: ",ProbaLeft,", rightInx: ",rightIndx,",ProbaRight: ",ProbaRight)
        if(verb>2)
          cat(", computational Bgd: ",timeGANMC," q: ",probaL$q)
        cat(". (Algorithm=",algo,")")
      }
      cat(". Memory currently used: ",sum(gc()[,2])," MB\n")
    }
    cat("\n")
  }

  ffSave=F
  j=0
  while(ProbaRight<alpha & (rightIndx-leftIndx)>=2){
    NextEval<-ceiling((rightIndx+leftIndx)/2)

    # if indMaxSort is less than 100 there's no need to use MCgrf
    if(rightIndx<sizeMaxGenz){
      if(type==">"){
        a<-threshold
        b<-Inf
      }else{
        a<--Inf
        b<-threshold
      }
      ProbaEval<-pmvnorm(lower = a,upper = b,mean = mu[1:NextEval],sigma = Sigma[1:NextEval,1:NextEval])
      varEval<-attr(ProbaEval,"error")
    }
    else{ # Otherwise we need to..
      precMC<-2.5e3
      qq<-c(min(15,leftIndx-1),min(25,leftIndx-1))
      if(!lightReturn)
        ffSave=T
      if((rightIndx-leftIndx)<=7){
        precMC<-5e3
        timeGANMC=25
        qq<-c(min(20,leftIndx-1),min(30,leftIndx-1))
        if((rightIndx-leftIndx)<=5){
          precMC<-1e4
          qq<-c(min(20,leftIndx-1),min(40,leftIndx-1))
        }
        if((rightIndx-leftIndx)<=3){
          precMC<-3e4
          timeGANMC=30
          qq<-c(min(25,leftIndx-1),min(60,leftIndx-1))
        }
      }

      if(type==">"){
        probaEval<-ProbaMin(cBdg=timeGANMC,q=qq,E=E[1:NextEval],threshold=threshold,mu=mu[1:NextEval],Sigma = Sigma[1:NextEval,1:NextEval],pn=NULL,lightReturn=!((!lightReturn)*ffSave),method=meth,verb=(verb-1),Algo=algo)
      }else{
        pn1E<-sortPn$x[1:NextEval]
        probaEval<-ProbaMax(cBdg=timeGANMC,q = qq,E = E[1:NextEval,],threshold = threshold,mu = mu[1:NextEval],Sigma = Sigma[1:NextEval,1:NextEval],pn = pn1E,lightReturn=!((!lightReturn)*ffSave),method=meth,verb=(verb-1),Algo=algo)
      }
      ProbaEval<-as.numeric(1-probaEval$probability)
      varEval<- probaEval$variance
    }


    if(ProbaEval>alpha){
      leftIndx<-NextEval
      ProbaLeft<-ProbaEval
      varLeft<-varEval
      rightIndx<-rightIndx
      ProbaRight<-ProbaRight
    }
    else{
      leftIndx<-leftIndx
      ProbaLeft<-ProbaLeft
      rightIndx<-NextEval
      ProbaRight<-ProbaEval
      varRight<-varEval
    }

    j=j+1
    if(verb){
      cat("\nIteration ",j," finished.\n")
      if(verb>1){
        if(rightIndx<sizeMaxGenz){
          cat("leftIndx:",leftIndx,",ProbaLeft: ",ProbaLeft,", rightInx: ",rightIndx,",ProbaRight: ",ProbaRight," (Genz)")
        }else{
          cat("leftIndx:",leftIndx,",ProbaLeft: ",ProbaLeft,", rightInx: ",rightIndx,",ProbaRight: ",ProbaRight)
          if(verb>2)
            cat(", timeGANMC: ",timeGANMC," q: ",qq)
          cat(" (Algo=",algo,")")
        }
        cat(". Memory currently used: ",sum(gc()[,2])," MB\n")
      }
      cat("\n")
    }

  }

  if(!lightReturn){
    resList<-list(set=(pn>sortPn$x[rightIndx]),lvs=sortPn$x[rightIndx],proba=ProbaLeft,vars=varLeft)#,lastMCit=finPpL)
  }else{
    resList<-list(set=(pn>sortPn$x[rightIndx]),lvs=sortPn$x[rightIndx])
  }
  if(verb){
    cat("Conservative estimate computation finished.")
  }
  rm(Sigma)
  if(verb){
    cat("(Memory currently used: ",sum(gc()[,2])," MB)\n")
  }
  return(resList)

}
