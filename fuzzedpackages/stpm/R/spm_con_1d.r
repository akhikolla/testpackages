#' Fitting a 1-D SPM model with constant parameters
#' 
#' @description This function implements a analytical solution to estimate the parameters in the continuous SPM model by assuming all the parameters are constants.
#' 
#' @param spm_data A dataset for the SPM model. See the STPM package for more details about the format.
#' @param a The initial value for the paramter \eqn{a}. The initial value will be predicted if not specified.
#' @param b The initial value for the paramter \eqn{b}. The initial value will be predicted if not specified.
#' @param q The initial value for the paramter \eqn{q}. The initial value will be predicted if not specified.
#' @param f The initial value for the paramter \eqn{f}. The initial value will be predicted if not specified.
#' @param f1 The initial value for the paramter \eqn{f_1}. The initial value will be predicted if not specified.
#' @param mu0 The initial value for the paramter \eqn{\mu_0} in the baseline hazard. The initial value will be predicted if not specified.
#' @param theta The initial value for the paramter \eqn{\theta} in the baseline hazard. The initial value will be predicted if not specified.
#' @param lower A vector of the lower bound of the parameters.
#' @param upper A vector of the upper bound of the parameters.
#' @param control A list of the control parameters for the optimization paramters.
#' @param global A logical variable indicating whether the MLSL (TRUE) or the L-BFGS (FALSE) algorithm is used for the optimization.
#' @param verbose A logical variable indicating whether initial information is printed.
#' @param ahessian A logical variable indicating whether the approximate (FALSE) or analytical (TRUE) Hessian is returned.
#' @return est The estimates of the parameters.
#' @return hessian The Hessian matrix of the estimates.
#' @return lik The minus log-likelihood.
#' @return con A number indicating the convergence. See the 'nloptr' package for more details.
#' @return message Extra message about the convergence. See the 'nloptr' package for more details.
#' @references He, L., Zhbannikov, I., Arbeev, K. G., Yashin, A. I., and Kulminski, A.M., 2017. Genetic stochastic process model for detecting pleiotropic and interaction effects with longitudinal data.
#' @export
#' @examples { 
#' library(stpm) 
#' dat <- simdata_cont(N=500)
#' colnames(dat) <- c("id", "xi", "t1", "t2", "y", "y.next")
#' res <- spm_con_1d(as.data.frame(dat), a=-0.05, b=2, q=1e-8, f=80, f1=90, mu0=1e-3, theta=0.08)
#'}
spm_con_1d <- function(spm_data, a = NA, b = NA, q = NA, f = NA, f1 = NA, mu0 = NA, theta = NA, lower = c(), upper = c(), control = list(xtol_rel = 1e-6), global = FALSE, verbose = TRUE, ahessian = FALSE)
{

  #a <- -0.1
  #b <- 1.2
  #q <- 0.001
  #f <- 25
  #f1 <- 26
  #mu0 <- 0.01
  #theta <- 0

  data <- spm_data
  data[which(is.na(data$y.next)),'y.next'] <- data[which(is.na(data$y.next)),'y']
  data$r0 <- 0
  aggdata <-aggregate(data[,c('id','t1','t2','xi')], by=list(data$id), FUN='max')
  tau <- aggdata$t2
  n_j <- as.data.frame(table(data$id))[,2]
  m0 <- data$y
  r0 <- data$r0
  yij <- data$y.next
  tij <- data$t2
  t0 <- data$t1
  delta <- aggdata$xi
  name_par <- c('a','b','q','f','f1','mu0','theta')

  if(is.na(f))
  {f <- mean(m0, na.rm=TRUE)}
  if(is.na(f1))
  {f1 <- mean(m0, na.rm=TRUE)}
  if(is.na(a))
  {a <- (-0.1)*abs(max(m0, na.rm=TRUE))/max(t0,na.rm=TRUE)}
  if(is.na(mu0))
  {mu0 <- sum(delta, na.rm=TRUE)/sum(tij-t0, na.rm=TRUE)}
  if(is.na(theta))
  {theta <- 0.00001}
  if(is.na(b))
  {b <- sqrt(abs(mean(yij-m0, na.rm=TRUE)))/abs(mean(tij-t0,na.rm=TRUE))}
  if(is.na(q))
  {q <- mu0/(0.25*f^2)}

  if(a>=0)
  {stop('The value for the argument \'a\' must be negative.')}
  if(b<0)
  {stop('The value for the argument \'b\' must be positive.')}
  if(q<0)
  {stop('The value for the argument \'q\' must be positive.')}
  #if(mu0<0)
  #{stop('The value for the argument \'mu0\' must be positive.')}

  param <- c(a,b,q,f,f1,mu0,theta)

  if(length(lower)==0)
  {lower <- c(a*5, b*0.1, q*0.01, f*0.5, f1*0.5, mu0*0.1, theta)}

  if(length(upper)==0)
  {upper <- c(a*0.05, b*5, q*10, f*2, f1*2, mu0*10, theta)}

  if(verbose == TRUE)
  {
	print('Initial values:')
	print(param)
	print('Lower bounds:')
	print(lower)
	print('Upper bounds:')
	print(upper)
  }

  # re <- lik_con_1d(param, m0,r0,tau,yij,delta,tij, n_j, t0)
  # re <- c(re, gr_con_1d(param, m0,r0,tau,yij,delta,tij, n_j, t0))
  if(global == FALSE)
  {
	re <- lbfgs( x0=param,fn=lik_con_1d,gr=gr_con_1d,m0=m0,r0=r0,tau=tau,yij=yij,delta=delta,tij=tij, n_j=n_j, t0=t0,lower=lower,upper=upper, control = control)
  }else{
	re <- mlsl( x0=param,fn=lik_con_1d,gr=gr_con_1d,m0=m0,r0=r0,tau=tau,yij=yij,delta=delta,tij=tij, n_j=n_j, t0=t0,lower=lower,upper=upper, control = control)
  }

  hessian <- optim(re$par, lik_con_1d, gr_con_1d, m0,r0,tau,yij,delta,tij, n_j, t0, method = "L-BFGS-B", lower = re$par, upper = re$par, hessian=TRUE)

  a_hes <- hessian$hessian
  if(ahessian == TRUE)
  {
	a_hes <- hes_loglik(re$par,m0,r0,tau,yij,delta,tij, n_j,t0)
  }
  
  colnames(a_hes) <- name_par
  rownames(a_hes) <- name_par

  #par_re <- matrix(re$par, 1, 7)
  #colnames(par_re) <- c('a','b','q','f','f1','mu0','theta')
  
  ##### Calculate p-values for estimates ####
  coef <- re$par
  names(coef) <- name_par
  fi <- solve(a_hes)   #### Fischer Information Matrix. Note we're not using the negative of the hessian here
  stderr <- sqrt(diag(fi))
  zscore <- coef/stderr
  pvalue <- 2*(1 - pnorm(abs(zscore)))
  results.par_re <- cbind(coef,stderr,zscore,pvalue)
  colnames(results.par_re) <- c("Coeff.", "Std. Err.", "z", "p value")
  

  list(est=results.par_re,hessian=a_hes, lik=re$value, con=re$convergence, message = re$message)
}
