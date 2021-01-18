setMethod("pcoordplot",
  signature(x = "IdtMclust"),
  function(x, title="Parallel Coordinate Plot", Seq=c("AllMidP_AllLogR","MidPLogR_VarbyVar"), model="BestModel", ...)
  {
    if (requireNamespace("GGally",quietly=TRUE)==FALSE) {
      stop("Required package GGally is not installed\n")
    }

    if (model=="BestModel") {
      G <- BestG(x)
      x_mean  <- mean(x)
    }  else {
      nomodmsg <- paste("There is no model named",model,"in these Idtmclust results\n") 
      if (substr(model,1,3)=="Hom") {
        if (!is.element(model,names(x@allres$RepresHom))) stop(nomodmsg) 
        G  <- x@allres$RepresHom[[model]]@nG
        x_mean  <- x@allres$RepresHom[[model]]@parameters$mean
      } else if (substr(model,1,3)=="Het") {
        if (!is.element(model,names(x@allres$RepresHet))) stop(nomodmsg)
        G  <- x@allres$RepresHet[[model]]@nG
        x_mean <- x@allres$RepresHet[[model]]@parameters$mean
      } else stop(nomodmsg) 
    }

    q <- x@NIVar   # Numbef or Interval-valued variables
    p <- 2*q       # Total number of MidPoints and LogRanges
    
    Seq <- match.arg(Seq)
    if (Seq == "AllMidP_AllLogR")  {
      DF <- as.data.frame(t(x_mean))
    }  else if (Seq == "MidPLogR_VarbyVar")  {
      DF <- as.data.frame(t(x_mean[rep(1:q,each=2)+rep(c(0,q),q),]))
    }
#    DF <- cbind(DF,Component=paste("CP",1:G,sep=""))
    DF <- cbind.data.frame(DF,Component=paste("CP",1:G,sep=""))
    plt <- ggparcoord(DF, columns = 1:p, groupColumn = 'Component') + ggtitle(title) +
      theme_minimal() + theme(plot.title=element_text(hjust=0.5)) + theme(axis.text.x=element_text(angle=90,hjust=1)) 
       
    plt <- plt + theme(axis.title = element_blank())

    plt <- plt + theme(axis.text.y = element_blank())
    plt <- plt + theme(panel.grid.major.y = element_blank())

    plt <- plt + theme(panel.grid.minor = element_blank())
    plt <- plt + theme(panel.grid.major.x = element_line(color = "#bbbbbb"))
    min_y <- min(plt$data$value)
    max_y <- max(plt$data$value)
    pad_y <- (max_y - min_y) * 0.1
    lab_x <- rep(1:p, times = 2) # 2 times, 1 for min 1 for max
    lab_y <- rep(c(min_y - pad_y, max_y + pad_y), each = p)

    z <- c(sapply(DF[,1:p], min), sapply(DF[, 1:p], max))
    absz <- abs(z)
    z_scinot <- which(absz<1e-3|absz>=1e4)
    z_3dig <- which(absz>=1e-3&absz<1.)
    z_2dig <- which(absz>=1.&absz<10.)
    z_1dig <- which(absz>=10.&absz<100.)
    z_0dig <- which(absz>=100.&absz<1e4)
    
    lab_z <- character(length(z))
    if (length(z_scinot)>0) {
      lab_z[z_scinot] <- formatC(z[z_scinot],digits=1,format="e")
    }
    if (length(z_3dig)>0) {
      lab_z[z_3dig] <- formatC(z[z_3dig],digits=3,format="f")
    }
    if (length(z_2dig)>0) {
      lab_z[z_2dig] <- formatC(z[z_2dig],digits=2,format="f")
    }
    if (length(z_1dig)>0) {
      lab_z[z_1dig] <- formatC(z[z_1dig],digits=1,format="f")
    }
    if (length(z_0dig)>0) {
      lab_z[z_0dig] <- formatC(z[z_0dig],digits=0,format="f")
    }

    plt <- plt + annotate("text", x = lab_x, y = lab_y, label = lab_z, size = 3)
    
    print(plt)
  }
)

setMethod("plotInfCrt",
  signature(object = "IdtMclust"),
  function(object, crt=object@SelCrit,legpos="bottomleft", nprnt=5, ...)
  {
    if (crt=="BIC") {
      IC_values <- object@BICs
    } else if (crt=="AIC") {
      IC_values <- object@AICs
    } else {
      stop("Wrong crt argument (it should equal either the string 'BIC' or the string 'AIC'.\n")
    }   
    if (length(IC_values)<2) stop("the argument to plotInfCrt does not include several (more than one) models to compare.\n")
    
    nn <- length(IC_values)
    modnames <- names(IC_values)
    Mxt <- substr(modnames,1,3)
    Mxts <- unique(Mxt)
    nMxts <- length(Mxts)
    nGdigts <- nchar(modnames) - 5
    G <- substr(names(IC_values),4,3+nGdigts)
    Gs <- unique(G)
    nGs <- length(Gs)
    CovCase <- substr(names(IC_values),3+nGdigts+1,nGdigts+6)
    CovCases <- unique(CovCase)
    nCovCases <- length(CovCases)
    MxbyCovC <- paste(rep(Mxts,each=nCovCases),rep(CovCases,nMxts),sep="") 
    nMxbyCovC <- nMxts*nCovCases
    if (nGs==1) {
      if (nMxts==2) IC_values <- IC_values[rep(1:nCovCases,each=2) + rep(c(0,nCovCases),nCovCases)]
      barplot(IC_values,ylab=crt,xlab="Model",space=10.,...)
      abline(h=0)
      return()  
    }  
    ICmat <- matrix(NA,nGs,nMxbyCovC,dimnames = list(Gs,MxbyCovC))
    linetype <- ifelse(substr(MxbyCovC,1,3)=="Hom",1,9)
    for (m in 1:nMxts) for (c in 1:nCovCases) 
#      ICmat[,(m-1)*nCovCases+c] <- IC_values[Mxt==Mxts[m] & CovCase==CovCases[c]] 
      ICmat[,(m-1)*nCovCases+c] <- -IC_values[Mxt==Mxts[m] & CovCase==CovCases[c]] 

#    matplot(as.integer(substr(Gs,2,1+unique(nGdigts))),ICmat,type="b",lty=1,ylab=crt,xlab="Number of Components",
    matplot(as.integer(substr(Gs,2,1+unique(nGdigts))),ICmat,type="b",lty=linetype,ylab=paste("minus",crt),xlab="Number of Components",
            pch=substr(MxbyCovC,5,5),col=1:nMxbyCovC,...)
#    if (max(ICmat[,1]) < max(ICmat[,nMxbyCovC])) legpos <- "topleft"
#    else  legpos <- "topright"
    legend(legpos,legend=MxbyCovC,lty=linetype,col=1:nMxbyCovC)

    bestcrt <- sort(IC_values,index.return=TRUE)
    cat("Best",crt,"values:\n")
    cat("          ",names(IC_values)[bestcrt$ix[1:nprnt]],"\n",sep="   ")
    cat(crt,"   ",bestcrt$x[1:nprnt],"\n",sep="   ")
    cat(paste(crt,"diff"),"     ",bestcrt$x[1:nprnt]-bestcrt$x[1],"\n",sep="   ")
  }
)  

setMethod("show",
  signature(object = "IdtMclust"),
  function(object)
  {
    cat("\'", class(object)[1], "\' model object:\n", sep = "")
    if (object@Hmcdt) {
      str1 <- "Homoscedastic"
    } else { 
      str1 <- "Heteroscedastic"
    }
    cat(" best model: ", str1, " setup with covariance configuration C",object@BestC," and ", object@BestG, " components\n",sep="") 
    invisible()
  }
)  

setMethod("summary",
  signature(object = "IdtMclust"),
  function(object, parameters = FALSE, classification = FALSE, model="BestModel", ...)
  {
    if (model=="BestModel") {
      G  <- object@BestG
      C  <- object@BestC
      mod <- object
    }  else {
      nomodmsg <- paste("There is no model named",model,"in these Idtmclust results\n") 
      if (substr(model,1,3)=="Hom") {
        if (!is.element(model,names(object@allres$RepresHom))) stop(nomodmsg) 
        mod <- object@allres$RepresHom[[model]]
      } else if (substr(model,1,3)=="Het") {
        if (!is.element(model,names(object@allres$RepresHet))) stop(nomodmsg)
        mod <- object@allres$RepresHet[[model]]
      } else stop(nomodmsg) 
      G  <- mod@nG
      C  <- mod@Conf
    }
    pro <- mod@parameters$pro
    if(is.null(pro)) pro <- 1
    names(pro) <- paste("CP",seq(G),sep="")
    mean <- mod@parameters$mean
    covariance <- mod@parameters$covariance
    title <- paste("Gaussian finite mixture model fitted by EM algorithm")
    if (mod@Hmcdt) {
      str1 <- "Homoscedastic"
    } else { 
      str1 <- "Heteroscedastic"
    }
    modelName <- paste(str1, " C",C,sep="") 

#    obj <- list(
#      title = title, modelName =modelName,  Hmcdt = mod@Hmcdt,     
#      NObs = mod@NObs, NIVar = mod@NIVar, G = G,  
#      loglik = mod@logLik, bic = mod@bic,
#      pro = pro, mean = mean, covariance = covariance,
#      classification = mod@classification,
#      printParameters = parameters, printClassification = classification
#    )
#    class(obj) <- "summaryIdtMclust"
#    return(obj)

   new("summaryIdtMclust",title=title,modelName=modelName,Hmcdt=mod@Hmcdt,     
      NObs=mod@NObs,NIVar=mod@NIVar,G=G,loglik=mod@logLik,bic=mod@bic,
      pro=pro,mean=mean,covariance=covariance,classification=mod@classification,
      printParameters=parameters,printClassification=classification)
  }
) 
 
#print.summaryIdtMclust <- function(x, ...)
#{
#   n <- x$NObs
#   p <- 2*x$NIVar
#   cat(rep("-", nchar(x$title)),"\n",sep="")
#   cat(x$title, "\n")
#   cat(rep("-", nchar(x$title)),"\n",sep="")
#   cat(x$modelName," model with ", x$G, ifelse(x$G > 1, " components\n", " component\n"),sep = "") 
#   tab <- data.frame("log-likelihood" = x$loglik, "NObs" = n, "BIC" = x$bic, row.names = "")
#   print(tab)
#   cat("\nClustering table:")
#   if (x$G==1) {
#     cat("\nCP1\n",n,"\n")
#   } else {
#     print(table(factor(x$classification, levels = paste("CP",seq_len(x$G),sep=""))))
#   }
#   if(x$printParameters) {
#     cat("\nMixing probabilities:\n")
#     print(x$pro)
#     cat("\nMeans:\n")
#     print(x$mean)
#     cat("\nStandard deviations:\n")
#     if (x$Hmcdt) {
#       print(sqrt(diag(x$covariance[,,1])))
#     } else {
#       stdv <- matrix(nrow=p,ncol=x$G,dimnames=list(rownames(x$mean),colnames(x$mean)))
#       for(g in 1:x$G) stdv[,g] <- sqrt(diag(x$covariance[,,g]))
#       print(stdv)
#     }
#     cat("\nCorrelations:\n")
#     if (x$Hmcdt) {
#       print(cov2cor(x$covariance[,,1]))
#     } else { 
#       for(g in 1:x$G) { 
#         cat("[,,CP", g, "]\n", sep = "")
#         print(cov2cor(x$covariance[,,g])) 
#       }
#     }
#   }
#   if(x$printClassification) {
#     cat("\nClassification:\n")
#     print(x$classification)
#   }
#   invisible(x)
#}

setMethod("show",
  signature(object = "summaryIdtMclust"),
  function(object)
  {
   n <- object@NObs
   p <- 2*object@NIVar
   cat(rep("-", nchar(object@title)),"\n",sep="")
   cat(object@title, "\n")
   cat(rep("-", nchar(object@title)),"\n",sep="")
   cat(object@modelName," model with ", object@G, ifelse(object@G > 1, " components\n", " component\n"),sep = "") 
   tab <- data.frame("log-likelihood" = object@loglik, "NObs" = n, "BIC" = object@bic, row.names = "")
   print(tab)
   cat("\nClustering table:")
   if (object@G==1) {
     cat("\nCP1\n",n,"\n")
   } else {
     print(table(factor(object@classification, levels = paste("CP",seq_len(object@G),sep=""))))
   }
   if(object@printParameters) {
     cat("\nMixing probabilities:\n")
     print(object@pro)
     cat("\nMeans:\n")
     print(object@mean)
     cat("\nStandard deviations:\n")
     if (object@Hmcdt) {
       print(sqrt(diag(object@covariance[,,1])))
     } else {
       stdv <- matrix(nrow=p,ncol=object@G,dimnames=list(rownames(object@mean),colnames(object@mean)))
       for(g in 1:object@G) stdv[,g] <- sqrt(diag(object@covariance[,,g]))
       print(stdv)
     }
     cat("\nCorrelations:\n")
     if (object@Hmcdt) {
       print(cov2cor(object@covariance[,,1]))
     } else { 
       for(g in 1:object@G) { 
         cat("[,,CP", g, "]\n", sep = "")
         print(cov2cor(object@covariance[,,g])) 
       }
     }
   }
   if(object@printClassification) {
     cat("\nClassification:\n")
     print(object@classification)
   }
   invisible(object)
  }
)

#Accessor methods

setMethod("parameters",signature(x = "IdtMclust"),
  function(x,model="BestModel") 
  { 
    if (model=="BestModel") return(x@parameters) 
    nomodmsg <- paste("There is no model named",model,"in these Idtmclust results\n") 
    if (substr(model,1,3)=="Hom") {
      if (!is.element(model,names(x@allres$RepresHom))) stop(nomodmsg) 
      return(x@allres$RepresHom[[model]]@parameters)
    } else if (substr(model,1,3)=="Het") {
      return(x@allres$RepresHet[[model]]@parameters)
    } else stop(nomodmsg) 
  }
)

setMethod("pro",signature(x = "IdtMclust"),
  function(x,model="BestModel") 
  { 
    if (model=="BestModel") return(x@parameters$pro) 
    nomodmsg <- paste("There is no model named",model,"in these Idtmclust results\n") 
    if (substr(model,1,3)=="Hom") {
      if (!is.element(model,names(x@allres$RepresHom))) stop(nomodmsg) 
      return(x@allres$RepresHom[[model]]@parameters$pro)
    } else if (substr(model,1,3)=="Het") {
      return(x@allres$RepresHet[[model]]@parameters$pro)
    } else stop(nomodmsg) 
  }
)

setMethod("mean",signature(x = "IdtMclust"),
  function(x,model="BestModel") 
  { 
    if (model=="BestModel") return(x@parameters$mean) 
    nomodmsg <- paste("There is no model named",model,"in these Idtmclust results\n") 
    if (substr(model,1,3)=="Hom") {
      if (!is.element(model,names(x@allres$RepresHom))) stop(nomodmsg) 
      return(x@allres$RepresHom[[model]]@parameters$mean)
    } else if (substr(model,1,3)=="Het") {
      return(x@allres$RepresHet[[model]]@parameters$mean)
    } else stop(nomodmsg) 
  }
)

setMethod("var",signature(x = "IdtMclust"),
  function(x) return(x@parameters$covariance)
)

setMethod("classification",signature(x = "IdtMclust"),
  function(x,model="BestModel") 
  { 
#    if (model=="BestModel") return(factor(x@parameters$classification)) 
    if (model=="BestModel") return(factor(x@classification)) 
    nomodmsg <- paste("There is no model named",model,"in these Idtmclust results\n") 
    if (substr(model,1,3)=="Hom") {
      if (!is.element(model,names(x@allres$RepresHom))) stop(nomodmsg) 
#      return(factor(x@allres$RepresHom[[model]]@parameters$classification))
      return(factor(x@allres$RepresHom[[model]]@classification))
    } else if (substr(model,1,3)=="Het") {
#      return(factor(x@allres$RepresHet[[model]]@parameters$classification))
      return(factor(x@allres$RepresHet[[model]]@lassification))
    } else stop(nomodmsg) 
  }
)

setMethod("SelCrit",signature(x = "IdtMclust"),function (x) x@SelCrit)
setMethod("Hmcdt",signature(x = "IdtMclust"),function (x) x@Hmcdt)
setMethod("BestG",signature(x = "IdtMclust"),function (x) x@BestG)
setMethod("BestC",signature(x = "IdtMclust"),function (x) x@BestC)
setMethod("PostProb",signature(x = "IdtMclust"),function(x) x@z)
setMethod("logLik",signature(object = "IdtMclust"),function(object,..) object@logLik)
setMethod("BIC",signature(object = "IdtMclust"),function(object,...) object@bic)
setMethod("AIC",signature(object = "IdtMclust"),function(object,...,k=2) object@aic)

setMethod("cor",signature(x ="IdtMclust"),
 function(x)
{
   if (x@Hmcdt==TRUE) {
     ncov <- 1
   } else {
     ncov <- x@BestG
   }
   covdim <- 2*x@NIVar
   covnames <- dimnames(x@parameters$covariance)[[1]]
   res <- array(dim=c(covdim,covdim,ncov),dimnames=list(covnames,covnames,NULL))
   for (g in 1:ncov) res[,,g] <- cov2cor(x@parameters$covariance[,,g])
   res
 }
) 



