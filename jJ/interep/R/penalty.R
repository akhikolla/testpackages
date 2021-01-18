#' This function gives the penalty functions
#' @importFrom stats gaussian
#' @param x matrix of covariates.
#' @param n the sample size.
#' @param p the number of predictors.
#' @param q the number of environment factors.
#' @param beta the coefficient vector.
#' @param lam1 the tuning parameter lambda1 for individual penalty.
#' @param pmethod the penalization method. "mixed" refers to MCP penalty to individual main effects and group MCP penalty to interactions; "individual" means MCP penalty to all effects.
#' @param p1 the number of gene factors.
#' @param lam2 the tuning parameter lambda2 for group penalty.
#' @return
#' \item{E}{the penalty function.}
#' @export

penalty <- function(x,n,p,q,beta,lam1,pmethod,p1,lam2){
  eps=0.000001
  if(pmethod=="individual"){
    beta.mcp=beta[(q+2):p]
    x.mcp=x[,(q+2):p]

    E.mcp=rep(0,(p-q-1))

    for (j in 1:(p-q-1)) {
      sub=j
      x.sub=x.mcp[,sub]
      beta0=beta.mcp[sub]
      kj=t(x.sub)%*%x.sub/n
      norm = sqrt(mean((x.sub*beta0)^2))
      #norm=as.numeric(sqrt(t(beta0)%*%kj%*%beta0))
      E.mcp[j]=dmcp(abs(as.vector(beta0)),lam1)/(abs(as.vector(norm))+eps)
      #E1=adiag(E1,E.groupmcp[,,j])
    }
    E=E.mcp
    E=c(rep(0,q+1),E)
    E=diag(E)
  }

  if(pmethod=="mixed"){
    beta.mcp=beta[1:(p1+q+1)]
    x.mcp=x[,1:(p1+q+1)]

    E.mcp=rep(0,(p1+q+1))

    for (j in 1:(p1+q+1)) {
      sub=j
      x.sub=x.mcp[,sub]
      beta0=beta.mcp[sub]
      kj=t(x.sub)%*%x.sub/n
      norm = sqrt(mean((x.sub*beta0)^2))
      E.mcp[j]=dmcp(abs(as.vector(beta0)),lam1)/(abs(as.vector(norm))+eps)
    }

    x.gmcp=x[,(p1+q+2):p]
    beta.gmcp=beta[c((p1+q+2):p)]
    for (j in 1:p1) {
      sub=((j-1)*q+1):(j*q)
      x.sub=x.gmcp[,sub]
      beta0=beta.gmcp[sub]
      norm = sqrt(mean((x.sub%*%beta0)^2))
      E.mcp=c(E.mcp,dmcp(abs(as.vector(beta0)),lam2)/(abs(as.vector(norm))+eps))
    }

    E1.mcp=diag(E.mcp)
    E1.mcp[,c(1:(1+q))]<-0
    E<-E1.mcp
  }

  return(E)
}
