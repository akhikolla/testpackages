## deflation function
Deflate.PCA <- function(x,v){
  ind <- which(v!=0)
  rx <- x
  if (length(ind) != 1){
    rx[,ind] <- x [,ind] - (x[,ind] %*% v[ind]) %*% t(v[ind])
  }else{
    rx[,ind] <- 0
  }
  return(rx)
}

### function to compute score outlying weights
# Input:
# x: data matrix;
# ws.od: logical variable to indicate whether an observation is an orthogonal outlier (ws.od=False) or not;
#        orthogonal outlierswill be removed from computing the score outlying weights.
# reference: P. Filzmoser, R. Maronna, M. Werner. Outlier identification in high dimensions, Computational Statistics and Data Analysis, 52, 1694-1711, 2008.
# "IdOUT" is rewritten from the fucntion "pcout" in "mvoutlier" package

IdOUT <- function(x,ws.od = ws.od,outbound=0.25){

  crit.M1 = 1/3
  crit.c1 = 2.5

  crit.M2 = 0.25
  crit.c2 = 0.99

  cs = 0.25

  x <- as.matrix(x)

  p1 <- dim(x)[2]

  x.pc <- as.matrix(x[ws.od,])
  xpc.sc <- scale(x,center =  apply(x.pc, 2, median), scale = apply(x.pc, 2, mad))

  ### the first weights
  ####### the kurtosis ##
  wp <- abs(apply(as.matrix(xpc.sc[ws.od,])^4, 2, mean) - 3)
  xpcw.sc <- xpc.sc %*% (diag(wp/sum(wp)))
  xpc.norm <- sqrt(apply(xpcw.sc^2, 1, sum))

  x.dist1 <- xpc.norm * sqrt(qchisq(0.5, p1))/median(xpc.norm[ws.od])
  M1 <- quantile(x.dist1[ws.od], crit.M1)
  const1 <- median(x.dist1[ws.od]) + crit.c1 * mad(x.dist1[ws.od])

  w1 <- (1 - ((x.dist1 - M1)/(const1 - M1))^2)^2
  w1[x.dist1 < M1] <- 1
  w1[x.dist1 > const1] <- 0

  ### the second weights
  xpc.norm <- sqrt(apply(xpc.sc^2, 1, sum))
  x.dist2 <- xpc.norm * sqrt(qchisq(0.5, p1))/median(xpc.norm[ws.od])
  M2 <- sqrt(qchisq(crit.M2, p1))
  const2 <- sqrt(qchisq(crit.c2, p1))
  w2 <- (1 - ((x.dist2 - M2)/(const2 - M2))^2)^2
  w2[x.dist2 < M2] <- 1
  w2[x.dist2 > const2] <- 0


  wfinal <- (w1 + cs) * (w2 + cs)/((1 + cs)^2)
  wfinal01 <- wfinal >= outbound

  return(list(w1 = w1, w2 = w2, wp = wp, wfinal = wfinal, wfinal01 = wfinal01, sd = x.dist2, wt.sd = x.dist1))
}



#################
# the following functions are originally written by Tom Reynkens in rospca packages
# Reference: Hubert M., Reynkens T., Schmitt E. & Verdonck T. Sparse PCA for High-Dimensional Data With Outliers, "Technometrics", 58, 424-434

matStand <- function (X, f_c = mean, f_s = sd, center= TRUE, scale = TRUE) {
  p <- ncol(X)
  if (center) {
    center <- apply(X,2,f_c)
  } else {
    center <- rep(0,p)
  }
  if (scale) {
    scale <- apply(X,2,f_s)
  } else {
    scale <- rep(1,p)
  }
  data <- sweep(X,2,center,"-")
  return(list(data=sweep(data,2,scale,"/"),c=center,s=scale))
}


unimcd <- function (x, h){
  nobs <- length(x)
  len <- nobs-h+1
  xorig <- x

  if (len==1) {
    tmcd <- mean(x)
    smcd <- sd(x)
    weights <- rep(1,nobs)

  } else {

    #Sort data and keep indices
    sorted <- sort(x,index.return=TRUE)
    x <- sorted$x
    I <- sorted$ix


    #####
    #Look for h-subset with smallest variance
    #We do this by looking at contiguous h-subsets which allows us
    #to compute the variances (and means) recursively

    #Sum of h-subsets
    ax <- numeric(len)
    ax[1] <- sum(x[1:h])
    for (samp in 2:len) {
      ax[samp] <- ax[samp-1]-x[samp-1]+x[samp+h-1]
    }

    #Sample variances of h-subsets multiplied by h-1
    ax2 <- ax^2/h
    sq <- numeric(len)
    sq[1] <- sum(x[1:h]^2)-ax2[1]
    for (samp in 2:len) {
      sq[samp] <- sq[samp-1]-x[samp-1]^2+x[samp+h-1]^2-ax2[samp]+ax2[samp-1];
    }

    #Subset(s) with minimal variance
    sqmin <- min(sq)
    ii <- which(sq==sqmin)
    ndup <- length(ii)
    optSets <- ax[ii]

    #initial mean is mean of h-subset
    initmean <- optSets[floor((ndup+1)/2)]/h
    #initial variance is variance of h-subset
    initcov <- sqmin/(h-1)

    #####
    # Consistency factor

    res <- (xorig-initmean)^2/initcov
    sortres <- sort(res)
    factor <- sortres[h]/qchisq(h/nobs,1)
    initcov <- factor*initcov

    #####
    #Weights and reweighting

    #raw squared robust distances
    res <- (xorig-initmean)^2/initcov
    quantile <- qchisq(0.975,1)

    #raw weights
    weights <- (res<=quantile)

    #reweighting procedure
    tmcd <- sum(xorig*weights)/sum(weights)
    smcd <- sqrt(sum((xorig-tmcd)^2*weights)/(sum(weights)-1))
  }

  return(list(tmcd=tmcd,smcd=smcd,weights=weights))

}



#Cutoff for the ODs, also returns ODs
coOD <- function (od, h = length(od)){

  #Problems can occur when the ODs are very small
  if (all(od<10^(-14))) {
    co.od <- 0; od <- rep(0,length(od))
  } else {
    mcd <- unimcd(od^(2/3),h=h)
    co.od <- (mcd$tmcd+mcd$smcd*qnorm(0.975))^(3/2)
  }
  return(list(od=od,co.od=co.od))
}



