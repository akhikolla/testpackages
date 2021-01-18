wCorrSim <- function(n, rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=FALSE, outstr="") {
  len <- max(c(length(n), length(rho), length(ML), length(fast), length(reset), length(usew)))
  vec <- c("n", "rho", "ML", "fast", "reset", "usew")
  for(i in 1:length(vec)) {
    var <- get(vec[i])
    if(length(var) != len) {
      if(length(var) != 1) {
        stop("length of ", sQuote(vec[i]), " must be 1 or the same as the longest vector passed to sim")
      } else {
        var <- rep(var,len)
      }
    }
    assign(vec[i],var)
  }
  everusew <- sum(usew)>0
  
  ns <- n
  cor0 <- rho
  df <- data.frame(n=n,rho=rho, ML=ML, fast=fast, reset=reset, usew=usew)
  

  df$spear <- df$speart <- NA
  df$Q <- df$M <- NA
  df$pear <- df$peart <- NA
  df$pc <- df$pct <- NA
  df$ps <- df$pst <- NA
  
  ii <- 1
  while(ii <= nrow(df)) {
    cori <- df$rho[ii]
    n <- df$n[ii]
    ML <- df$ML[ii]
    fast <- df$fast[ii]
    reset <- df$reset[ii]
    usew <- df$usew[ii]

    if(interactive()) {
      cat(outstr,"n=",n,"cori=",cori,"pct=",100*ii/nrow(df),"\n")
      cat("  fast=",fast,"ml=",ML,"reset=",reset,"\n")
    }
    if(reset) {
      n <- ifelse(everusew, 10*df$n[ii], df$n[ii])

      cr <- cori
      x <- y <- w <- M <- Q <- c()
      while(length(w) < df$n[ii]) {
        xp <- rnorm(n)
        yp <- sqrt(1-cr^2)*rnorm(n) + cr*xp
        if(everusew) {
          wp <- (xp-yp)^2+1
          pr <- 1/wp
          pr <- df$n[ii] * pr/(sum(pr) * 100)
          wp <- 1/pr
          samp <- (1:n)[runif(n)<pr]
          x <- c(x,xp[samp])
          y <- c(y,yp[samp])
          w <- c(w,wp[samp])
        } else {
          x <- xp
          y <- yp
          w <- rep(1/n, n)
        }
      }
      M <- 1  
      Q <- 1
      nm <- sample(2:5,1)
      nq <- sample(2:5,1)
      x <- x[1:df$n[ii]]
      y <- y[1:df$n[ii]]
      w <- w[1:df$n[ii]]
      iter <- 1
      while(  ((length(unique(M)) < 2) | (length(unique(Q)) < 2)) & (iter < 100)) {
        iter <- iter + 1
        tm <- sort(rnorm(nm))
        tq <- sort(rnorm(nq))
        theta1 <- c(NA,-Inf,tq,Inf)
        theta2 <- c(NA,-Inf,tm,Inf)

        Q <- rep(NA,n)
        for(i in 2:length(theta1)) {
          Q <- ifelse(x>theta1[i], i, Q)
        }
        Q <- Q - 1
        Q <- as.numeric(as.factor(Q))

        M <- rep(NA,n)
        for(i in 2:length(theta2)) {
          M <- ifelse(y>theta2[i], i, M)
        }
        M <- M - 1
        M <- as.numeric(as.factor(M))
      }
      if(iter >=99) {
        cat("could not get multiple bins\n")
        cat("x <- c(",paste(x,collapse=","),")\n")
        cat("y <- c(",paste(y,collapse=","),")\n")
        cat("M <- c(",paste(M,collapse=","),")\n")
        cat("Q <- c(",paste(Q,collapse=","),")\n")
      }
      df$M[ii] <- length(unique(M))
      df$Q[ii] <- length(unique(Q))
    } else {
      df$cor[ii] <- df$cor[ii-1]
      df$M[ii] <- length(unique(M))
      df$Q[ii] <- length(unique(Q))
    }
    
    if(usew) {
      wu <- w
    } else {
      wu <- rep(1,length(x))
    }
    
    st0 <- system.time(fcorp <- weightedCorr(x,y, method="Pearson", weights=wu, fast=fast, ML=ML))
    df$peart[ii] <- sum(st0[1:2])
    df$pear[ii] <- fcorp

    st0 <- system.time(fcorp <- weightedCorr(x,y, method="Spearman", weights=wu, fast=fast, ML=ML))
    df$speart[ii] <- sum(st0[1:2])
    df$spear[ii] <- cor(x,y)
    
    st0 <- system.time(fcorp <- weightedCorr(x,M, method="Polyserial", weights=wu, fast=fast, ML=ML))
    df$pst[ii] <- sum(st0[1:2])
    df$ps[ii] <- fcorp
    
    st0 <- system.time(fcorp <- weightedCorr(M, Q, method="Polychoric", weights=wu, fast=fast, ML=ML))
    df$pct[ii] <- sum(st0[1:2])
    df$pc[ii] <- fcorp
    
    ii <- ii + 1
  }
  dfout <- data.frame(n=rep(df$n,4),
                      rho=rep(df$rho,4),
                      ML=rep(df$ML,4),
                      usew=rep(df$usew,4),
                      fast=rep(df$fast,4),
                      est=c(df$pear, df$spear, df$ps, df$pc),
                      t=c(df$peart, df$speart, df$pst, df$pct),
                      type=rep(c("Pearson", "Spearman", "Polyserial", "Polychoric"),each=nrow(df)))
  dfout
}

spearmanSim <- function(n, rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=TRUE, outstr="") {
  len <- max(c(length(n), length(rho), length(ML), length(fast), length(reset), length(usew)))
  vec <- c("n", "rho", "ML", "fast", "reset", "usew", "Spearman", "N")
  N <- n
  Spearman <- rho
  for(i in 1:length(vec)) {
    var <- get(vec[i])
    if(length(var) != len) {
      if(length(var) != 1) {
        stop("length of ", sQuote(vec[i]), " must be 1 or the same as the longest vector passed to sim")
      } else {
        var <- rep(var,len)
      }
    }
    assign(vec[i],var)
  }
  everusew <- sum(usew)>0
  
  ns <- n
  cor0 <- rho
  df <- data.frame(n=n,rho=rho, ML=ML, fast=fast, reset=reset, usew=usew, Spearman = rho)
  
  
  df$spear <- df$speart <- NA
  df$pear <- df$peart <- NA
  
  ii <- 1
  while(ii <= nrow(df)) {
    cori <- df$rho[ii]
    n <- df$n[ii]
    ML <- df$ML[ii]
    fast <- df$fast[ii]
    reset <- df$reset[ii]
    usew <- df$usew[ii]
    
    if(interactive()) {
      cat(outstr,"n=",n,"cori=",cori,"pct=",100*ii/nrow(df),"\n")
      cat("  fast=",fast,"ml=",ML,"reset=",reset,"\n")
    }
    if(reset) {
      n <- ifelse(everusew, 50*df$n[ii], df$n[ii])
      
      cr <- cori
      x <- y <- w <- M <- Q <- c()
      while(length(w) < df$n[ii]) {
        xp <- rnorm(n)
        yp <- sqrt(1-cr^2)*rnorm(n) + cr*xp
        df$Spearman[ii] <- cor(xp, yp, method="spearman")
        if(everusew) {
          wp <- (xp-yp)^2+1
          pr <- 1/wp
          pr <- pr*df$n[ii]/(sum(pr))
          wp <- 1/pr
          #samp <- sample(1:n, size=df$n[ii], replace=FALSE, prob=pr)
          samp <- (1:n)[runif(n)<pr]
          #print(samp)
          df$N[ii] <- length(samp)
          x <- c(x,xp[samp])
          y <- c(y,yp[samp])
          w <- c(w,wp[samp])
        } else {
          x <- xp
          y <- yp
          w <- rep(1/n, n)
        }
      }
      x <- x[1:df$n[ii]]
      y <- y[1:df$n[ii]]
      w <- w[1:df$n[ii]]
    } else {
      df$cor[ii] <- df$cor[ii-1]
    }
    
    if(usew) {
      wu <- w
    } else {
      wu <- rep(1,length(x))
    }
    
    
    
    st0 <- system.time(fcorp <- weightedCorr(x,y, method="Spearman", weights=wu, fast=fast, ML=ML))
    df$speart[ii] <- sum(st0[1:2])
    df$spear[ii] <- fcorp
    
    
    ii <- ii + 1
  }
  dfout <- data.frame(n=rep(df$n,1),
                      rho=rep(df$rho,1),
                      ML=rep(df$ML,1),
                      usew=rep(df$usew,1),
                      fast=rep(df$fast,1),
                      est=c(df$spear),
                      t=c(df$speart),
                      SpearmanOrig = df$Spearman,
                      N = df$N,
                      type=rep(c("Spearman"),each=nrow(df)))
  dfout
}




