#'@title Variable Coefficients for Zeta
#'@aliases summary_precision
#'@name summary_precision
#'@description A function that uses the gamma values of the posterior distribution of the model and calculates the estimates for each zeta covariate.
#'@usage summary_precision(x,prob=0.95)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param prob a probability containing the credibility index for the HPD interval for the coefficients of the covariates.
#'@return A list containing the estimates for the covariables of zeta, this list contains the following items:
#'\describe{
#'\item{table}{a table with the means, medians, standard deviations and the Highest Posterior Density (HPD) Interval,}
#'\item{coeff}{a vector containing the estimated coefficients for the variables.}}
#'@seealso \code{\link{summary_mean}},\code{\link{values}},\code{\link{summary.bayesbr}}
summary_precision = function(x,prob=0.95){
  gamma = x$info$samples$gamma
  zetas = x$info$samples$zeta
  names_w = x$info$names$names_w
  table = NULL
  coeff = numeric()
  if(!is.null(gamma)){
    table = c()
    tam = length(gamma)
    for (i in 1:tam) {
      if(i==1 & tam==1){
        gammas = gamma$'gammas'
      }
      else{
        aux = paste0('gammas[',i,']')
        gammas = gamma[[aux]]
      }
      gamma_mcmc = as.mcmc( c(gammas) )
      hpd = HPDinterval(gamma_mcmc, prob=prob)
      mean_t = round(mean(gammas),5)
      coeff = c(coeff,mean_t)
      median_t = round(median(gammas),5)
      sd_t = round(sd(gammas),5)
      vec = c(mean_t,median_t,sd_t,round(hpd[1:2],5))
      table = rbind(table,vec)
    }
    table = matrix(table,ncol = 5)
    colnames(table) = c("Mean","Median", "Std. Dev.","HPD_inf","HPD_sup")
    rownames(table) = names_w
    names(coeff) = names_w
  }
  else{
    zeta = zetas$zeta
    zeta_mcmc = as.mcmc( c(zeta) )
    hpd = HPDinterval(zeta_mcmc, prob=prob)
    mean_t = round(mean(zeta),5)
    coeff = c(coeff,mean_t)
    median_t = round(median(zeta),5)
    sd_t = round(sd(zeta),5)
    vec = c(mean_t,median_t,sd_t,round(hpd[1:2],5))
    table = rbind(table,vec)
    table = matrix(table,ncol = 5)
    colnames(table) = c("Mean","Median", "Std. Dev.","HPD_inf","HPD_sup")
    rownames(table) = c("(phi)")
    names(coeff) = c("(phi)")
  }
  list = list(table = table,gammas = coeff)
  return(list)
}
