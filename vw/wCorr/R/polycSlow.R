# based losely on Olsson, Ulf (1979), "Maximum Likelihood Estimation of the Polychoric Correlation Coefficient", Psychometrica, 44(4), 443-460.
#' @importFrom mnormt biv.nt.prob
#' @importFrom minqa bobyqa
#' @importFrom stats qnorm
#' @importFrom stats optimize
#' @importFrom stats cor
polycSlow <- function(x,y,w,ML=FALSE) {

  lnl <- function(xytab, cc, rc, corr) {
    cc <- c(-Inf, cc, Inf)
    rc <- c(-Inf, rc, Inf)
    pm <- sapply(1:(length(cc)-1), function(c) {
            sapply(1:(length(rc)-1), function(r) {
              biv.nt.prob(df=Inf,
                          lower=c(cc[c], rc[r]),
                          upper=c(cc[c+1], rc[r+1]),
                          mean=c(0,0),
                          S=matrix(c(1,corr,corr,1), nrow=2, ncol=2, byrow=TRUE))
              #pmvnorm(lower=c(cc[c], rc[r]),
              #        upper=c(cc[c+1], rc[r+1]),
              #        mean=c(0,0),
              #        corr=matrix(c(1,corr,corr,1), nrow=2, ncol=2, byrow=TRUE))
          })
    })
    suppressWarnings(lpm <- log(pm))
    #lpm[is.nan(lpm)] <- 0
    lpm[(is.nan(lpm)) | (!is.finite(lpm))] <- log(.Machine$double.xmin)
    sum(xytab * lpm)
  }

  optf_all <- function(par, xytab) {
    c1 <- ncol(xytab)-1
    c2 <- c1 + nrow(xytab)-1
    -1 * lnl(xytab, cc=fscale_cuts(par[1:c1]), rc=fscale_cuts(par[(c1+1):c2]), corr=fscale_corr(par[length(par)] ))
  }

  optf_corr <- function(par, xytab, theta1, theta2) {
    c1 <- ncol(xytab)-1
    c2 <- c1 + nrow(xytab)-1
    -1 * lnl(xytab, cc=fscale_cuts(theta2), rc=fscale_cuts(theta1), corr=fscale_corr(par))
  }

  fscale_cuts <- function(par) {
    cumsum(c(par[1],exp(par[-1])))
  }

  fscale_corr <- function(par) {
    tanh(par)
  }

  weightedTable <- function(x,y,w=rep(1,length(x))) {
    tab <- table(x,y)
    for(i in 1:nrow(tab)) {
      for(j in 1:ncol(tab)) {
        tab[i,j] <- sum(w[ x==dimnames(tab)[[1]][i] & y == dimnames(tab)[[2]][j] ])
      }
    }
    tab
  }

  imapTheta <- function(theta0) {
    c(theta0[1], log(theta0[-1]-theta0[-length(theta0)]))
  }
  
  xytab <- weightedTable(x,y,w)

  # first check for perfect correlations which throw the optimizer for a loop because of the infinite bounds of the mapped correlation
  i <- 1
  j <- 1
  foundConcord <- FALSE
  foundDiscord <- FALSE
  while(j<ncol(xytab)) {
    if(i<nrow(xytab) & j < ncol(xytab)) {
      if(xytab[i,j]>0 & sum(xytab[(i+1):nrow(xytab), (j+1):ncol(xytab)]) > 0) {
        foundConcord <- TRUE
        break
      }
    }
    if(i>1 & j > 1) {
      if(xytab[i,j]>0 & sum(xytab[1:(i-1), 1:(j-1)]) > 0) {
        foundConcord <- TRUE
        break
      }
    }
    # incriment
    i <- i + 1
    if(i>nrow(xytab)) { 
      i <- 1
      j <- j + 1
    }
  }

  i <- 1
  j <- 1
  while(j<ncol(xytab)) {
    if(i>1 & j < ncol(xytab)) {
      if(xytab[i,j]>0 & sum(xytab[1:(i-1), (j+1):ncol(xytab)]) > 0) {
        foundDiscord <- TRUE
        break
      }
    }
    if(i<nrow(xytab) & j > 1) {
      if(xytab[i,j]>0 & sum(xytab[(i+1):nrow(xytab), 1:(j-1)]) > 0) {
        foundDiscord <- TRUE
        break
      }
    }
    # incriment
    i <- i + 1
    if(i>nrow(xytab)) { 
      i <- 1
      j <- j + 1
    }
  }
  if(!foundDiscord){
    #print(xytab)
    return(1)
  }
  if(!foundConcord) {
    #print(xytab)
    return(-1)
  }
 
  #GKgamma <- rcorr.cens(x,y,outx=T)["Dxy"]
  #if( GKgamma %in%  c(-1,1)) {
  #  return(unname(GKgamma))
  #}

  #op <- optim(par=c(log(1:(ncol(xytab)-1)), log(1:(nrow(xytab)-1)),cor(x,y)), optf_all, xytab=xytab, control=list(fnscale=-1), method="BFGS")
  #fscale_corr(op$par[length(op$par)])
  ux <- sort(unique(x))
  cut1 <- imapTheta( sapply(ux[-length(ux)],function(z) qnorm(sum(w[x<=z])/sum(w)) ))
  uy <- sort(unique(y))
  cut2 <- imapTheta( sapply(uy[-length(uy)],function(z) qnorm(sum(w[y<=z])/sum(w)) ))
  
  cor0 <- atanh(cor(as.numeric(x),as.numeric(y)))
  #bob <- bobyqa(c(cut1,cut2,cor0), fn=optf_all, xytab=xytab)
  if(ML) {
    bob <- bobyqa(c(cut1,cut2,cor0), fn=optf_all, xytab=xytab)
    return(fscale_corr(bob$par[length(bob$par)]))
  } else {
    opt <- optimize(optf_corr, interval=cor0+c(-3,3), xytab=xytab, theta1=cut1,theta2=cut2)
    return(  fscale_corr(opt$minimum))
  }
  # should return above
}
