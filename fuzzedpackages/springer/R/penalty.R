#' This function provides the penalty functions.  Users can choose one of the three penalties: sparse group MCP, group MCP and MCP.
#' @importFrom stats gaussian
#' @param x the matrix of predictors, consisting of the clinical covariates, environmental factors, genetic factors and gene-environment interactions.
#' @param n the sample size.
#' @param t the number of clinical covariates.
#' @param p the number of predictors, which consists of the clinical covariates, environmental factors, genetic factors and gene-environment interactions.
#' @param q the number of environment factors.
#' @param beta the coefficient vector.
#' @param lam1 the tuning parameter \eqn{\lambda_1} for individual-level penalty.
#' @param structure Three choices are available for structured variable selection. "bilevel" for sparse-group selection on both group-level and individual-level. "group" for selection on group-level only. "individual" for selection on individual-level only.
#' @param p1 the number of genetic factors.
#' @param lam2 the tuning parameter \eqn{\lambda_2} for group-level penalty.
#' @return
#' \item{H}{the penalty function.}
#' @details
#' When {structure="bilevel"}, sparse group MCP is adopted and variable selection for longitudinal data including both genetic main effects and gene-environment interactions will be conducted on both individual and group levels (bi-level selection):
#' \itemize{
#' \item \strong{Group-level selection:} If the \eqn{v}th genetic factor has any effect at all (associated with the response or not) can be determined by whether \eqn{||\eta_{v}||_{2}=0}.
#' \item \strong{Individual-level selection:} whether the \eqn{v}th genetic variant has main effect, G\eqn{\times}E interaction or both can be determined by the nonzero componet.
#' }
#' If {structure="group"}, group MCP will be used and only group-level selection will be conducted on \eqn{||\eta_{v}||_{2}}; if {structure="individual"}, MCP will be adopted and only individual-level selection will be conducted on each \eqn{\eta_{vu}}, (\eqn{u=1,\ldots,q}).
#'
#' The minimax concave penalty (MCP) is adopted as the baseline penalty function in the springer package. Methods based on other popular choices, such as SCAD and LASSO, will be examined in the future.
#' @export

penalty <- function(x,n,t,p,q,beta,lam1,structure,p1,lam2){
  eps=0.000001
  if(structure=="individual"){
    beta.mcp=beta[(t+q+2):p]
    x.mcp=x[,(t+q+2):p]

    E.mcp=rep(0,(p-q-t-1))

    for (j in 1:(p-q-t-1)) {
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
  }

  if(structure=="group"){
    x.gmcp=x[,(t+q+2):p]
    beta.gmcp=beta[(t+q+2):p]
    E.gmcp=rep(0,(p-q-t-1))

    for (j in 1:p1) {
      sub=((j-1)*(q+1)+1):(j*(q+1))
      x.sub=x.gmcp[,sub]
      beta0=beta.gmcp[sub]
      norm = sqrt(mean((x.sub%*%beta0)^2))
      # kj=t(x.sub)%*%x.sub/n
      # norm=as.numeric(sqrt(t(beta0)%*%kj%*%beta0))
      E.gmcp[sub]=dmcp(abs(as.vector(beta0)),lam2)/(abs(as.vector(norm))+eps)
      # E.groupmcp[,,j]=diag(q_mcp(abs(as.vector(beta0)),lam2)/(abs(as.vector(norm))+eps))
      # E1=adiag(E1,E.groupmcp[,,j])
    }
    E=E.gmcp
  }

  if(structure=="bilevel"){
    beta.mcp=beta[(t+q+2):p]
    x.mcp=x[,(t+q+2):p]

    E.mcp=rep(0,(p-q-t-1))

    for (j in 1:(p-q-t-1)) {
      sub=j
      x.sub=x.mcp[,sub]
      beta0=beta.mcp[sub]
      kj=t(x.sub)%*%x.sub/n
      norm = sqrt(mean((x.sub*beta0)^2))
      #norm=as.numeric(sqrt(t(beta0)%*%kj%*%beta0))
      E.mcp[j]=dmcp(abs(as.vector(beta0)),lam1)/(abs(as.vector(norm))+eps)
      #if(j%%(q+1)==1) E.mcp[j]=0
    }

    x.gmcp=x[,(t+q+2):p]
    beta.gmcp=beta[(t+q+2):p]
    E.gmcp=rep(0,(p-q-t-1))

    for (j in 1:p1) {
      sub=((j-1)*(q+1)+1):(j*(q+1))
      x.sub=x.gmcp[,sub]
      beta0=beta.gmcp[sub]
      norm = sqrt(mean((x.sub%*%beta0)^2))
      # kj=t(x.sub)%*%x.sub/n
      # norm=as.numeric(sqrt(t(beta0)%*%kj%*%beta0))
      E.gmcp[sub]=dmcp(abs(as.vector(beta0)),lam2)/(abs(as.vector(norm))+eps)
      # E.groupmcp[,,j]=diag(q_mcp(abs(as.vector(beta0)),lam2)/(abs(as.vector(norm))+eps))
      # E1=adiag(E1,E.groupmcp[,,j])
    }
    E=E.gmcp+E.mcp
  }

  E=c(rep(0,t+q+1),E)
  E=diag(E)
  H=E
  return(H)
}
