# based losely on Olsson, Ulf (1979), "Maximum Likelihood Estimation of the Polychoric Correlation Coefficient", Psychometrica, 44(4), 443-460.
#' @importFrom mnormt biv.nt.prob
#' @importFrom minqa bobyqa
#' @importFrom stats qnorm
#' @importFrom stats optimize
#' @importFrom stats cor
polycFast <- function(x,y,w,ML=FALSE) {

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
      })
    })
    lnlFast(xytab, pm)
  }

  optf_all <- function(par, xytab) {
    c1 <- ncol(xytab)-1
    c2 <- c1 + nrow(xytab)-1
    -1 * lnl(xytab, cc=fscale_cutsFast(par[1:c1]), rc=fscale_cutsFast(par[(c1+1):c2]), corr=fscale_corr(par[length(par)] ))
  }

  optf_corr <- function(par, xytab, theta1, theta2) {
    c1 <- ncol(xytab)-1
    c2 <- c1 + nrow(xytab)-1
    -1 * lnl(xytab, cc=fscale_cutsFast(theta2), rc=fscale_cutsFast(theta1), corr=fscale_corr(par))
  }


  fscale_corr <- function(par) {
    tanh(par)
  }
  xytab <- tableFast(x,y,w)

  #GKgamma <- rcorr.cens(x,y,outx=T)["Dxy"]

 # if (!(GKgamma %in%  c(-1,1)))
 #   if (discord %in% c(-1,1))
 #   {
 #     print(paste("Discord = ", discord(xytab)))
 #     print(paste("GK = ", unname(GKgamma)))
 #   }

  # if( GKgamma %in%  c(-1,1)) {
  #  return(unname(GKgamma))
  # }


  temp <- discord(xytab)

  if(temp==-1 | temp == 1)
     return(temp)

  ux <- sort(unique(x))
  
  cut1 <- imapThetaFast( sapply(ux[-length(ux)],function(z) qnorm(sum(w[x<=z])/sum(w)) ))
  uy <- sort(unique(y))
  cut2 <- imapThetaFast( sapply(uy[-length(uy)],function(z) qnorm(sum(w[y<=z])/sum(w)) ))
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
