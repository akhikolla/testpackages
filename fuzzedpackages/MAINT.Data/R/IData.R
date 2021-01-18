setMethod("summary",
  signature(object = "IData"),
  function (object) 
  {
    Rngsumar <- summary(exp(object@LogR))
    dimnames(Rngsumar)[[2]] <- paste(object@VarNames,".Range",sep="")

#cat("summary.IData 1 -- got here\n")
#cat("str(object@MidP) =\n") ; print(str(object@MidP))
#cat("str(summary(object@MidP)) =\n") ; print(str(summary(object@MidP)))
#cat("str(Rngsumar =\n") ; print(str(Rngsumar))
#cat("str(object@LogR) =\n") ; print(str(object@LogR))
#cat("str(summary(object@LogR)) =\n") ; print(str(summary(object@LogR)))

#cat("methods(print) =\n") ; print(methods(print))
#cat("getS3method(print,summaryIData,optional = TRUE) =\n") ; getS3method("print","summaryIData",optional = TRUE)

#cat("summary.IData 2 -- got here\n")


#    new("summaryIData",MidPsumar=summary(object@MidP),Rngsumar=Rngsumar,LogRsumar=summary(object@LogR))
    show(new("summaryIData",MidPsumar=summary(object@MidP),Rngsumar=Rngsumar,LogRsumar=summary(object@LogR)))
  }
)

setMethod("show",
  signature(object = "IData"),
  function (object) 
  {
    printrow <- function(Bnds,NIVar) 
    { 
      cat(Bnds[1],"  ")
      for (j in 2:(NIVar+1))
        cat("[",format(Bnds[j],width=8,digits=5,justify="centre"),", ",
          format(Bnds[NIVar+j],width=8,digits=5,justify="centre"),"]  ",sep="") 
      cat("\n")
    }

    HalfRange <- exp(object@LogR)/2
    LB <- object@MidP - HalfRange
    UB <- object@MidP + HalfRange
    lobsname <- max(nchar(object@ObsNames))
    flength <- max(nchar(object@VarNames),nchar(format(LB[1,1],width=8,digits=5))+nchar(format(UB[1,1],width=8,digits=5))) + 6 
    for (j in 1:object@NIVar) {
      if (j>1) {
        nspaces <- flength-nchar(object@VarNames[j])
      } else {
        nspaces <- lobsname+ceiling(flength/2)-ceiling(nchar(object@VarNames[1])/2)+3
      }  
      cat(rep(" ",nspaces),object@VarNames[j],sep="" )
    }  
    cat("\n") 
    apply(cbind.data.frame(format(object@ObsNames,width=lobsname),LB,UB),1,printrow,NIVar=object@NIVar)
    invisible(object)
  }
)

setMethod("nrow",signature(x = "IData"),function(x) x@NObs)
setMethod("ncol",signature(x = "IData"),function(x) x@NIVar)
setMethod("dim",signature(x = "IData"),function(x) c(nrow(x),ncol(x)))
setMethod("rownames",signature(x = "IData"),function(x) x@ObsNames)
setMethod("colnames",signature(x = "IData"),function(x) x@VarNames)
setMethod("names",signature(x = "IData"),function(x) x@VarNames)
setMethod("MidPoints",signature(Idt = "IData"),function(Idt) Idt@MidP)
setMethod("LogRanges",signature(Idt = "IData"),function(Idt) Idt@LogR)
setMethod("Ranges",signature(Idt = "IData"),function(Idt) exp(Idt@LogR))

setMethod("NbMicroUnits",
  signature(x = "IData"),
  function(x) {
    if (length(x@NbMicroUnits)==0) return(NULL)
    x@NbMicroUnits
  }
)


setMethod("head",
  signature(x = "IData"),
  function (x,n=min(nrow(x),6L)) 
  {
    
    if (n>0) {
#      x[1:n,] 
      x[1:n,1:x@NIVar] 
    } else {
#      x[-1:n,]
      x[-1:n,1:x@NIVar]
    }
  }
)

setMethod("tail",
  signature(x = "IData"),
  function (x,n=min(nrow(x),6L)) 
  {
    if (n>0) {
#      x[(x@NObs-n+1):x@NObs,] 
      x[(x@NObs-n+1):x@NObs,1:x@NIVar] 
    } else {
#      x[(-x@NObs-n-1):-x@NObs,]
      x[(-x@NObs-n-1):-x@NObs,1:x@NIVar]
    }
  }
)

setMethod("plot",
  signature(x = "IData",y = "IData"),
  function(x, y, type=c("crosses","rectangles"), append=FALSE,  ...)
  {
    if (x@NIVar > 1) stop("Currently IData method plot can plot only one integer variable on the horizontal axis\n")
    if (y@NIVar > 1) stop("Currently IData method plot can plot only one integer variable on the vertical axis\n")
    if (x@NObs != y@NObs) stop("Arguments x and y have a different number of elements\n")

    type <- match.arg(type)

    dotarguments <- match.call(expand.dots=FALSE)$...

    if (is.null(dotarguments$main)) {
      dotarguments$main <- paste(y@VarNames,"vs.",x@VarNames)
    }
    if (is.null(dotarguments$xlab)) {
      dotarguments$xlab <- x@VarNames
    }
    if (is.null(dotarguments$ylab)) {
      dotarguments$ylab <- y@VarNames
    }


    HlfRngx <- exp(x@LogR[,1])/2
    Lbx <- x@MidP[,1] - HlfRngx
    Ubx <- x@MidP[,1] + HlfRngx
    HlfRngy <- exp(y@LogR[,1])/2
    Lby <- y@MidP[,1] - HlfRngy
    Uby <- y@MidP[,1] + HlfRngy

    if (is.null(dotarguments$xlim)) {
      minvx <- min(Lbx)
      maxvx <-max(Ubx)
      dotarguments$xlim <- c(0.95*minvx,1.05*maxvx)
    }
    if (is.null(dotarguments$ylim)) {
      minvy <- min(Lby)
      maxvy <-max(Uby)
      dotarguments$ylim <- c(0.95*minvy,1.05*maxvy)
    }

    if (!append) {
#      do.call("plot",c(list(x=0.,y=0.,type="n"),dotarguments))
      do.call("plot.default",c(list(x=0.,y=0.,type="n"),dotarguments))
    }

    if (type=="crosses") {
      do.call("segments",c(list(x0=Lbx,y0=y@MidP[,1],x1=Ubx,y1=y@MidP[,1]),dotarguments))
      do.call("segments",c(list(x0=x@MidP[,1],y0=Lby,x1=x@MidP[,1],y1=Uby),dotarguments))
    } else if (type=="rectangles") {
      for (i in 1:x@NObs) rect(Lbx[i],Lby[i],Ubx[i],Uby[i],lty=1)
    }
  }
)

setMethod("plot",
  signature(x = "IData",y = "missing"),
  function(x, casen=NULL, layout=c("vertical","horizontal"), append=FALSE,  ...)
  {
    if (x@NIVar > 1) {
      if (x@NIVar==2) {
        plot.default(x[,1],x[,2],...)
        return()
      } else {
        stop("Currently method plot can only plot at most tow interval variables simultaneously\n")
      }
    }

    layout <- match.arg(layout)

    dotarguments <- match.call(expand.dots=FALSE)$...
    if (is.null(dotarguments$main)) {
      dotarguments$main <- x@VarNames
    }
    if (is.null(dotarguments$ylab)) {
#      dotarguments$ylab <- x@VarNames
      if (layout=="vertical") dotarguments$ylab <- x@VarNames
      else dotarguments$ylab <- "Case numbers"
    }
    if (is.null(dotarguments$xlab)) {
#        dotarguments$xlab <- "Case numbers"
        if (layout=="vertical") dotarguments$xlab <- "Case numbers"
        else dotarguments$xlab <- x@VarNames
    }

    HlfRngy <- exp(x@LogR[,1])/2
    Lby <- x@MidP[,1] - HlfRngy
    Uby <- x@MidP[,1] + HlfRngy
    xcord <- 1:x@NObs
    if (!append) xcord <- 1:x@NObs
    else xcord <- 1:x@NObs + runif(x@NObs,-0.5,0.5)

    if (is.null(dotarguments$ylim) && layout=="vertical") {
      minvy <- min(Lby)
      maxvy <-max(Uby)
      dotarguments$ylim <- c(0.95*minvy,1.05*maxvy)
    }
    if (is.null(dotarguments$xlim) && layout=="horizontal") {
      minvy <- min(Lby)
      maxvy <-max(Uby)
      dotarguments$xlim <- c(0.95*minvy,1.05*maxvy)
    }

    if (is.null(casen)) {
      casen <- c(0,x@NObs+1)
    } else {
      casen <- factor(eval(casen))
    }
    if (layout=="vertical") {
      if (!append) {
        do.call("plot.default",c(list(x=casen,y=rep(0.,length(casen)),type="n"),dotarguments))
      }
      do.call("segments",c(list(x0=xcord,y0=Lby,x1=xcord,y1=Uby),dotarguments))
    } else  {
      if (!append) {
        do.call("plot.default",c(list(x=rep(0.,length(casen)),y=casen,type="n"),dotarguments))
      }
      do.call("segments",c(list(y0=xcord,x0=Lby,y1=xcord,x1=Uby),dotarguments))
    }
  }
)

setMethod("show",
  signature(object = "summaryIData"),
  function(object)
  {

#cat("show(summaryIData) -- got here")
#cat("str(object@MidPsumar) =",str(object@MidPsumar),"\n") 
#cat("str(object@MidPsumar) =",str(object@Rngsumar),"\n") 
#cat("str(object@MidPsumar) =",str(object@LogRsumar),"\n") 

    cat("Mid-Points summary:\n") ; print(object@MidPsumar)
    cat("Ranges summary:\n") ; print(object@Rngsumar)
    cat("Log-Ranges summary:\n") ; print(object@LogRsumar)
    invisible(object)
  }
)

#print.summaryIData <- function(x,...) 
#{
#  cat("Mid-Points summary:\n") ; print(x$MidPsumar)
#  cat("Ranges summary:\n") ; print(x$Rngsumar)
#  cat("Log-Ranges summary:\n") ; print(x$LogRsumar)
#  invisible(x)
#}

setMethod("print", signature(x="summaryIData"), function(x,...) invisible(x) )

IData <- function(Data,
  Seq=c("LbUb_VarbyVar","MidPLogR_VarbyVar","AllLb_AllUb","AllMidP_AllLogR"),VarNames=NULL,ObsNames=row.names(Data))
{
  if ( !is.data.frame(Data) && !is.matrix(Data) ) stop("First argument of IData must be a data frame or a matrix\n")
  p <- ncol(Data)  # Total number of Interval variable bounds
  q <- p/2	 # Number of Interval variables
  if (floor(q) != q) stop("Number of columns of Data ( =",p,") must be an even number\n")
  Seq <- match.arg(Seq)
  if (  (Seq == "LbUb_VarbyVar") || (Seq == "AllLb_AllUb") )
  {
    if (Seq == "LbUb_VarbyVar") { Lbnd <- Data[,2*(0:(q-1))+1] ; Ubnd <- Data[,2*(1:q)] }
    if (Seq == "AllLb_AllUb")   { Lbnd <- Data[,1:q] ; Ubnd <- Data[,(q+1):p] }
    MidP <- (Lbnd+Ubnd)/2
    LogR <- log(Ubnd-Lbnd)
    if (any(is.na(LogR))) stop("Invalid data")
  } else {
    if (Seq == "MidPLogR_VarbyVar") { MidP <- Data[,2*(0:(q-1))+1] ; LogR <- Data[,2*(1:q)] }
    if (Seq == "AllMidP_AllLogR")   { MidP <- Data[,1:q] ; LogR <- Data[,(q+1):p] }
  }
  if (is.null(VarNames)) VarNames <- paste("I",1:q,sep="")
  if (!is.data.frame(MidP)) MidP <- as.data.frame(MidP)
  if (!is.data.frame(LogR)) LogR <- as.data.frame(LogR)
  names(MidP) <- paste(VarNames,".MidP",sep="")
  names(LogR) <- paste(VarNames,".LogR",sep="")
  rownames(MidP) <- rownames(LogR) <- ObsNames

  new("IData",MidP=MidP,LogR=LogR,ObsNames=ObsNames,VarNames=VarNames,NObs=nrow(MidP),NIVar=q,NbMicroUnits=integer(0))
}

# Standard operators for IData objects

# Indexing and assignement

setMethod("[",
  signature(x='IData',i='ANY',j='ANY'),
  function(x,i,j,...,drop=TRUE)
  {
    if (missing(i)) i <- 1:nrow(x)  
    if (missing(j)) j <- 1:ncol(x)  
    if (class(i)[1]=="character") i <- sapply(i,function(ri) which(ri==x@ObsNames))
    if (class(j)[1]=="character") j <- sapply(j,function(ci) which(ci==x@VarNames))
    IData(cbind.data.frame(x@MidP[i,j,drop=FALSE],x@LogR[i,j,drop=FALSE]),
      Seq="AllMidP_AllLogR",VarNames=x@VarNames[j],ObsNames=x@ObsNames[i])
  }
)

setMethod("[<-",
  signature(x='IData',i='ANY',j='ANY'),
  function(x,i,j,value)
  {
    if (!is(value,"IData")) stop("Argument value is not an IData object\n")
    if (missing(i)) i <- 1:nrow(x)  
    if (missing(j)) j <- 1:ncol(x)  
    if (class(i)[1]=="character") i <- sapply(i,function(ri) which(ri==x@ObsNames))
    if (class(j)[1]=="character") j <- sapply(j,function(ci) which(ci==x@VarNames))
    x@MidP[i,j] <- value@MidP
    x@LogR[i,j] <- value@LogR
    x
  }
)

# Comparison

setMethod("==",
  signature(e1 = "IData",e2 = "IData"),
  function(e1,e2)
  {
    CompIvalue <- function(Ival) Ival[1,1] == Ival[1,2] && Ival[2,1] == Ival[2,2]

    if (!is(e2,"IData")) 
      stop("Trying to compare an IData object with an object of a diferent type\n")
    if ( e1@NObs != e2@NObs || e1@NIVar != e2@NIVar )
      stop("== only defined for equally-sized IData objects\n")
    TmpArray <- array(dim=c(e1@NObs,e1@NIVar,2,2))
    for (j in 1:e1@NIVar)  {
#      TmpArray[,j,,1] <- cbind(e1@MidP[,j],e1@LogR[,j])
#      TmpArray[,j,,2] <- cbind(e2@MidP[,j],e2@LogR[,j])
      TmpArray[,j,,1] <- cbind.data.frame(e1@MidP[,j],e1@LogR[,j])
      TmpArray[,j,,2] <- cbind.data.frame(e2@MidP[,j],e2@LogR[,j])
    }

    apply(TmpArray,c(1,2),CompIvalue)
  }
)

setMethod("!=",
  signature(e1 = "IData",e2 = "IData"),
  function(e1,e2)  
  {
    if (!is(e2,"IData"))
      stop("Trying to compare an IData object with an object of a diferent type\n")
    if ( e1@NObs != e2@NObs || e1@NIVar != e2@NIVar )
      stop("!= only defined for equally-sized IData objects\n")
    !(e1==e2)
  }
)

setMethod(rbind,
  signature(x = "IData"),
  function(x, y, ...)
  {
    if (class(x)[1]!="IData") stop("Argument x is not an object of class IData\n")
    if (class(y)[1]!="IData") stop("Argument y is not an object of class IData\n")
    if (x@NIVar != y@NIVar) stop("Arguments x and y have a different number of interval-valued variables\n")
#    dataDF <- rbind(cbind(x@MidP,x@LogR),cbind(y@MidP,y@LogR))
    dataDF <- rbind.data.frame(cbind.data.frame(x@MidP,x@LogR),cbind.data.frame(y@MidP,y@LogR))
    if (x@NIVar==1) {
      ONames <- c(x@ObsNames,y@ObsNames)
    } else {
      ONames <- rownames(dataDF)
    }
    par <- match.call()
    lpar <- length(par)
    if (lpar>4) for (i in 3:(lpar-2))  {
      newIDtObj <- eval(par[[i]])
      if (i==3) auxtxt <- "rd " else auxtxt <- "-th "
      if (class(newIDtObj)[1]!="IData") stop("The ",i,auxtxt,"argument is not an object of class IData\n")
      if (x@NIVar != newIDtObj@NIVar) stop("First and ",i,auxtxt,"argument have a different number of interval-valued variables\n")
#      dataDF <- rbind(dataDF,cbind(newIDtObj@MidP,newIDtObj@LogR))
      dataDF <- rbind.data.frame(dataDF,cbind.data.frame(newIDtObj@MidP,newIDtObj@LogR))
      if (x@NIVar==1) {
        ONames <- c(ONames,newIDtObj@ObsNames)
      } else {
        ONames <- rownames(dataDF)
      }
    }  
    IData(dataDF,Seq="AllMidP_AllLogR",ObsNames=ONames,VarNames=x@VarNames)
  }
)

setMethod(cbind,
  signature(x = "IData"),
  function(x, y, ...)
  {
    if (class(x)[1]!="IData") stop("Argument x is not an object of class IData\n")
    if (class(y)[1]!="IData") stop("Argument y is not an object of class IData\n")
    if (x@NObs != y@NObs) stop("Arguments x and y have a different number of rows\n")
#    dataDF <- cbind(x@MidP,y@MidP,x@LogR,y@LogR)
    dataDF <- cbind.data.frame(x@MidP,y@MidP,x@LogR,y@LogR)
    if (x@NIVar==1 && y@NIVar==1) rownames(dataDF) <- x@ObsNames
    VNames <- c(x@VarNames,y@VarNames)
    par <- match.call()
    lpar <- length(par)
    curnIvar <- x@NIVar + y@NIVar
    if (lpar>4) for (i in 3:(lpar-2))  {
      newIDtObj <- eval(par[[i]])
      if (i==3) auxtxt <- "rd " else auxtxt <- "-th "
      if (class(newIDtObj)[1]!="IData") stop("The ",i,auxtxt,"argument is not an object of class IData\n")
      if (x@NObs != newIDtObj@NObs) stop("First and ",i,auxtxt,"argument have a different number of rows\n")
#      dataDF <- cbind(dataDF[,1:curnIvar],newIDtObj@MidP,dataDF[,(curnIvar+1):(2*curnIvar)],newIDtObj@LogR)
      dataDF <- cbind.data.frame(dataDF[,1:curnIvar],newIDtObj@MidP,dataDF[,(curnIvar+1):(2*curnIvar)],newIDtObj@LogR)
      VNames <- c(VNames,newIDtObj@VarNames)
      curnIvar <- curnIvar + newIDtObj@NIVar
    }  
    IData(dataDF,Seq="AllMidP_AllLogR",VarNames=VNames)
  }
)
