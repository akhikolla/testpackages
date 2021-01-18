#' Bivariate Alternating Recurrent Event Response and Covariate Data Simulation
#'
#' @description This function simulates a series of alternating recurrent events based on the simulation setting in Lee, Huang, Xu, Luo (2018).
#'
#' @param nsize Sample size which refers to the number of subjects in the data set where each subject could have multiple episodes of events.
#' @param beta1 True coefficients for Type I gap times in the accelerated failure time model (AFT).
#' @param beta2 True coefficients for Type II gap times in the accelerated failure time model (AFT).
#' @param tau_c Maximum support of censoring time. It can take values as follows:
#' \itemize{
#' \item \verb{tau_c=63} (default): corresponds to a 15\% censoring rate for each scenario in Tables 1 and 2 of Lee, Huang, Xu, Luo (2018).
#' \item \verb{tau_c=30}: corresponds to a 30\% censoring rate for each scenario in Tables 1 and 2 of Lee, Huang, Xu, Luo (2018).
#' }
#'
#' @param set Simulation setting based on scenarios outlined in Tables 1 and 2 in Lee, Huang, Xu, Luo (2018). Choose 1.1 (default) for scenario 1 with \eqn{\rho=1} in the covariance matrix of the frailty vector, 1.2 for scenario 1 with \eqn{\rho=0.5}, 1.3 for scenario 1 with \eqn{\rho=0} and 2.0 for scenario 2.
#'
#' @return Data frame with the alternating recurrent event data and one continuous and one binary covariate.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom utils tail
#'
#'
#' @examples
#' library(BivRec)
#' set.seed(1234)
#' sim_data <- simBivRec(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5),
#'             tau_c=63, set=1.1)
#' head(sim_data)
#'
#' @references
#'  Lee CH, Huang CY, Xu G, Luo X. (2018). Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine, 37: 996-1008.
#' \url{https://doi.org/10.1002/sim.7563}

#' @export

##-----data generation
simBivRec <- function(nsize, beta1, beta2, tau_c, set) {

  if (missing(tau_c)) {tau_c <- 63}
  if (missing(set)) {set <- 1.1}
  sg2 <- 0.5

  id=1:nsize

  ##generate covariates (A1,A2)
  A1=rbinom(nsize,1,0.5)
  A2=runif(nsize)
  A=cbind(A1,A2)

  ##generate gamma
  if (set==1.1) {
    Sig=matrix(c(sg2,sqrt(sg2)*sqrt(sg2)*1,sqrt(sg2)*sqrt(sg2)*1,sg2),2,2)
    gamma=mvrnorm(nsize,c(1,1),Sig)
    if (nsize > 1) {
      gamma1=gamma[,1]
      gamma2=gamma[,2]} else {
        gamma1=gamma[1]
        gamma2=gamma[2]
      }
  } else {
    if (set==1.2) {
      Sig=matrix(c(sg2,sqrt(sg2)*sqrt(sg2)*0.5,sqrt(sg2)*sqrt(sg2)*0.5,sg2),2,2)
      gamma=MASS::mvrnorm(nsize,c(1,1),Sig)
      if (nsize > 1) {
        gamma1=gamma[,1]
        gamma2=gamma[,2]} else {
          gamma1=gamma[1]
          gamma2=gamma[2]
        }
    } else {
      if (set==1.3) {
        gamma1=rnorm(nsize,1,sqrt(sg2))
        gamma2=rnorm(nsize,1,sqrt(sg2))
      } else {
        if (set==2.0) {
          gamma1=rnorm(nsize,1,sqrt(sg2))
          gamma2=rgamma(nsize,shape=1/sg2,rate=1/sg2)
        } else {
          stop("Invalid value for set. Must be 1.1, 1.2, 1.3, or 2.0.")
        }
      }
    }
  }

  dat=NULL
  deltas=NULL

  if (tau_c != 63) {
    if (tau_c != 30) {
     stop("Invalid value for tau_c. Must be 63 or 30.")
    }
  }

  for (i in id) {
    x.tmp=exp(gamma1[i]+A[i,]%*%beta1)
    y.tmp=exp(gamma2[i]+A[i,]%*%beta2)

    ci=runif(1,0,tau_c)

    sum.z.tmp=0
    dat.tmp=NULL
    while(sum.z.tmp<ci) {
      #generate alternating gap times until C_i
      xij=x.tmp*exp(rnorm(1,0,sqrt(0.1)))
      yij=y.tmp*exp(rnorm(1,0,sqrt(0.1)))

      dat.tmp=rbind(dat.tmp,c(xij,yij,sum.z.tmp))
      sum.z.tmp=sum.z.tmp+(xij+yij)
    }
    epi=nrow(dat.tmp)
    #set censoring indicator and last censored gap time pair
    if (epi==1) {cen=cbind(1,0)} else {
      cen=cbind(rep(1,epi),c(rep(1,epi-1),0))}
    if (dat.tmp[epi,1]>(ci-dat.tmp[epi,3])) {
      dat.tmp[epi,1]=ci-dat.tmp[epi,3]
      dat.tmp[epi,2]=0
      cen[epi,1]=0
    } else {
      dat.tmp[epi,2]=ci-dat.tmp[epi,3]-dat.tmp[epi,1]
    }
    subdat=cbind(i,1:epi,dat.tmp[,1],dat.tmp[,2],ci, cen,matrix(A[i,],epi,2,byrow=TRUE))
    dat=rbind(dat,subdat)
    deltas=rbind(deltas, cen)
  }
  dat=as.data.frame(dat)
  colnames(dat)=c("id","epi","xij","yij","ci","d1", "d2","a1","a2")

  ##return data, censoring rate, average number of gap time pairs
  # cenrate=sum(table(dat$id)==1)/nsize
  # mbar=mean(table(dat$id))
  # maxm=max(table(dat$id))
  return(data=dat)
}
