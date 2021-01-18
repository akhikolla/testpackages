setMethod("show",
  signature(object = "IdtOutl"),
  function(object) {
    print(object@outliers)
    invisible(object)
  }
)  

setMethod("plot",
  signature(x = "IdtOutl",y = "missing"),
  function(x,scale=c("linear","log"),RefDist=getRefDist(x),eta=geteta(x),multiCmpCor=getmultiCmpCor(x),...)
  {
    scale <- match.arg(scale)
    MD2 <- getMahaD2(x)
    n <- x@NObs
    p <- x@p
    h <- x@h

    oneminuseta <- 1-eta
    if (multiCmpCor=="always") {
      oneminusalpha <- oneminuseta^(1/n)
    } else if (multiCmpCor=="never") {
      oneminusalpha <- oneminuseta
    }

    if (RefDist=="ChiSq")  {
      delta <- qchisq(oneminusalpha,p)
    } else if (RefDist=="HardRockeAdjF")  {
      delta <- qHardRoqF(oneminusalpha,n,p,h)
    } else if (RefDist=="HardRockeAsF")  {
      delta <- qHardRoqF(oneminusalpha,n,p,h,adj=FALSE)
    } else if (RefDist=="CerioliBetaF")  {
      delta1 <- ((h-1)^2/h) * qbeta(oneminusalpha,p/2,(h-p-1)/2)
      delta2 <- (((h+1)*(h-1)*p)/(h*(h-p))) * qf(oneminusalpha,p,h-p)
    }

   if (RefDist!="CerioliBetaF")  {
     if (scale=="linear") {
#       plot(MD2,main="Robust Mahalanobis Distances",xlab="",xaxt="n",...)
       plot.default(MD2,main="Robust Mahalanobis Distances",xlab="",xaxt="n",...)
     } else if (scale=="log") {
#       plot(MD2,main="Robust Mahalanobis Distances (log scale)",xlab="",xaxt="n",log="y",...)
       plot.default(MD2,main="Robust Mahalanobis Distances (log scale)",xlab="",xaxt="n",log="y",...)
     }
     axis(1,1:n,labels=names(MD2),las=2,...)
     abline(h=delta)
   } else {
      RewSet <- x@boolRewind
      if (scale=="linear") {
#        plot(MD2,main="Robust Mahalanobis Distances",xlab="",xaxt="n",type="n",...)
        plot.default(MD2,main="Robust Mahalanobis Distances",xlab="",xaxt="n",type="n",...)
      } else if (scale=="log") {
#        plot(MD2,main="Robust Mahalanobis Distances (log scale)",xlab="",xaxt="n",type="n",log="y",...)
        plot.default(MD2,main="Robust Mahalanobis Distances (log scale)",xlab="",xaxt="n",type="n",log="y",...)
        axis(1,1:n,labels=names(MD2),las=2,...)
      }
      axis(1,1:n,labels=names(MD2),las=2,...)      
      RewSetInd <- which(RewSet==TRUE)
      UnRewSetInd <- which(RewSet==FALSE)
      points(x=RewSetInd,y=MD2[RewSet],pch=19,col="blue")
      points(x=UnRewSetInd,y=MD2[!RewSet],pch=19,col="red")
      abline(h=delta1,col="blue")
      abline(h=delta2,col="red")
   }
  }
)


setMethod("getMahaD2",signature(IdtOtl = "IdtOutl"),function(IdtOtl) IdtOtl@MD2) 
setMethod("geteta", signature(IdtOtl = "IdtOutl"), function(IdtOtl) IdtOtl@eta)
setMethod("getRefDist", signature(IdtOtl ="IdtOutl"), function(IdtOtl) IdtOtl@RefDist)
setMethod("getmultiCmpCor", signature(IdtOtl ="IdtOutl"), function(IdtOtl) IdtOtl@multiCmpCor)

