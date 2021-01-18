#############################################
#####  Plot coefficients' paths         #####
#############################################
plot.l1mstateCoef = function(x, trans = NULL,...){
  p = x$numcovs
  lambda = x$fit[,1]
  nlambda = length(x$fit[,1])
  Q = x$numtrans
  if(max(trans)>Q){
    warning("Invalid transition value, it should not be larger the total transitions.")
    return()
  }else{
    if(length(trans) == 0){
      warning("Plot all the coefficient paths of all transitions")
      rescov = list()
      for(q in 1:Q){
        rescov[[q]] = matrix(rep(0,p*nlambda), nrow = p)
        for(i in 1:nlambda){
          rescov[[q]][,i] = x$aBetaO[[i]][q,]
        }
        plotCoef(as.matrix(rescov[[q]]), lambda = lambda)
      }
    }else{
      rescov = list()
      for(q in trans){
        rescov[[q]] = matrix(rep(0,p*nlambda), nrow = p)
        for(i in 1:nlambda){
          rescov[[q]][,i] = x$aBetaO[[i]][q,]
        }
        plotCoef(as.matrix(rescov[[q]]), lambda = lambda)
      }
    }
  }
}

plotCoef=function(beta,lambda,label=TRUE,xlab="log(Lambda)",ylab="Coefficients", title=NULL){
  which <- nonzeroCoef(beta)
  nwhich<-length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={
           warning("No plot produced since all coefficients zero")
           return()
         },
         "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  beta<-as.matrix(beta[which,,drop=FALSE])
  index<-log(lambda)
  col_set <- rainbow_hcl(nrow(beta))
  matplot(index,t(beta),
          type = "l", lty = "solid", lwd = 3, ylab = "Coefficients",
          xlab = "log(Lambda)", main = title, col = col_set,
          cex.axis = 1.25, font = 2, cex.lab = 1.5, col.lab = '#993333', font.lab=2)
  
  if(label){
    nnz=length(which)
    xpos=min(index)
    pos=3
    xpos=rep(xpos,nnz)
    ypos=beta[,ncol(beta)]
    text(xpos,ypos-.05,paste(which),cex=.8,pos=pos)
  }
}

nonzeroCoef = function (beta, bystep = FALSE) 
{
  ### bystep = FALSE means which variables were ever nonzero
  ### bystep = TRUE means which variables are nonzero for each step
  nr<-nrow(beta)
  if (nr == 1) {#degenerate case
    if (bystep) 
      apply(beta, 2, function(x) if (abs(x) > 0) 
        1
        else NULL)
    else {
      if (any(abs(beta) > 0)) 
        1
      else NULL
    }
  }
  else {
    beta<-abs(beta)>0 # this is sparse
    which<-seq(nr)
    ones<-rep(1,ncol(beta))
    nz<-as.vector((beta%*%ones)>0)
    which<-which[nz]
    if (bystep) {
      if(length(which)>0){
        beta<-as.matrix(beta[which,,drop=FALSE])
        nzel <- function(x, which) if (any(x)) 
          which[x]
        else NULL
        which<-apply(beta, 2, nzel, which)
        if(!is.list(which))which<-data.frame(which)# apply can return a matrix!!
        which
      }
      else{
        dn<-dimnames(beta)[[2]]
        which<-vector("list",length(dn))
        names(which)<-dn
        which
      }
    }
    else which
  }
}