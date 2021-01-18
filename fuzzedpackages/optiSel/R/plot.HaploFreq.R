
"plot.HaploFreq"<-function(x, ID=1, hap=1, refBreed=NULL, Chr=NULL, show.maxFreq=FALSE, ...){
  Freq <- x
  ismax <- FALSE
  if("freq" %in% names(Freq)){
    ismax      <- length(attributes(Freq)$refBreeds)>1
    refBreed   <- paste(attributes(Freq)$refBreeds, collapse = ", ")
    Freq       <- freqlist(Freq)
    names(Freq)<- refBreed
  }
  if(is.null(refBreed)){refBreed<-names(Freq)[1]}
  if(!(refBreed%in%names(Freq))){stop("This reference breed is not available.\n")}
  if(is.null(Chr)){Chr <- unique(attributes(Freq)$map$Chr)}
  map       <- attributes(Freq)$map
  Indiv     <- attributes(Freq)$Indiv
  thisBreed <- attributes(Freq)$thisBreed

  for(b in names(Freq)){Freq[[b]]<-Freq[[b]][map$Chr%in%Chr,]}
  map <- map[map$Chr%in%Chr, ]
  
  x    <- cumsum(tapply(map$Mb,map$Chr,max)[as.character(unique(map$Chr))])
  xlab <- (x+c(0,x[-length(x)]))/2
  x    <- setNames(c(0, x[-length(x)]), names(x))
  map$x<- map$Mb+x[as.character(map$Chr)]

  if(is.numeric(ID)){
    Nr <- ID
    k  <- 2*(Nr-1)+hap
    ID <- Indiv[2*Nr]
  }else{
    Nr <- (which.max(Indiv==ID)+1)/2
    k  <- 2*(Nr-1)+hap
  }
  Info    <- paste("of the ", thisBreed, " with ID ",ID," (Indiv No ",Nr,", Hap ", k, ")",sep="")
  maxFreq <- do.call(pmax,Freq)

  plot(map$x, 0*map$x, type="n",ylim=c(0,1), main="Haplotype Frequency in Reference Breeds",xlab="Marker Position in bp", xaxt="n", bty="n", ylab="Frequency", ...)
  y   <- apply(maxFreq,1,max)
  col <- rep(c("grey","darkgrey"),length(Chr))

  if(show.maxFreq){
    for(i in 1:length(Chr)){
      x <- map$x[map$Chr==Chr[i]]
      polygon(c(min(x),x,max(x)), c(0, y[map$Chr==Chr[i]],0), col=col[i], border=col[i])
    }
  }
  
  for(i in 1:length(Chr)){
    x <- map$x[map$Chr==Chr[i]]
    polygon(c(min(x),x,max(x)), c(0, Freq[[refBreed]][map$Chr==Chr[i],k],0), col="red", border="red")
  }

  lines(c(0,0),c(0,1))
  for(i in 1:length(Chr)){
    x <- map$x[map$Chr==Chr[i]]
    lines(max(x)*c(1,1),c(0,1))
    lines(x, maxFreq[map$Chr==Chr[i],k])
  }
  lines(c(0,max(map$x)),c(1,1))
  mtext(Info, cex=0.8)
  axis(1,at=xlab,labels=paste("Chr", names(xlab)),lty=0)
  lines(c(0,max(map$x)),c(0,0))
  if(length(Freq)==1){
    if(ismax){legend(-0.02*max(map$x[map$Chr==Chr[1]]),1.0,c(paste("max in", names(Freq))),lty=1, col=c("black"),cex=0.8, bty="n", adj=c(0,0.0))}
    if(!ismax){legend(-0.02*max(map$x[map$Chr==Chr[1]]),1.0,c(paste(   "in", names(Freq))),lty=1, col=c("black"),cex=0.8, bty="n", adj=c(0,0.0))}
  }else{
    legend(-0.02*max(map$x[map$Chr==Chr[1]]),1.0,c(paste("max in", paste(names(Freq),collapse=", ")), paste("in", refBreed)),lty=1, col=c("black", "red"),cex=0.8, bty="n", adj=c(0,0.0))
  }
}
