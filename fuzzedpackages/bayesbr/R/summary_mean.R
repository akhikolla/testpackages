#'@title Variable Coefficients for Theta
#'@aliases summary_mean
#'@name summary_mean
#'@description A function that uses the beta values of the posterior distribution of the model and calculates the estimates for each theta covariate.
#'@usage summary_mean(x,prob=0.95)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param prob a probability containing the credibility index for the HPD interval for the coefficients of the covariates.
#'@return A list containing the estimates for the covariables of theta, this list contains the following items:
#'\describe{
#'\item{table}{a table with the means, medians, standard deviations and the Highest Posterior Density (HPD) Interval,}
#'\item{coeff}{a vector containing the estimated coefficients for the variables.}}
#'@seealso \code{\link{summary_precision}},\code{\link{values}},\code{\link{summary.bayesbr}}
summary_mean = function(x,prob=0.95){
  beta = x$info$samples$beta
  names_x = x$info$names$names_x
  warmup = x$info$warmup
  iter = x$info$iter
  table = NULL
  coeff = numeric()
  if(!is.null(beta)){
    table = c()
    tam = length(beta)
    for (i in 1:tam) {
      if(i==1 && tam==1){
        betas = beta$'betas'
      }
      else{
        aux = paste0('betas[',i,']')
        betas = beta[[aux]]
      }
      mean_t = round(mean(betas),5)
      coeff = c(coeff,mean_t)
      median_t = round(median(betas),5)
      sd_t = round(sd(betas),5)
      beta_mcmc = as.mcmc( c(betas) )
      hpd = HPDinterval(beta_mcmc, prob=prob)
      vec = c(mean_t,median_t,sd_t,round(hpd[1:2],5))
      table = rbind(table,vec)
    }
    colnames(table) = c("Mean","Median", "Std. Dev.","HPD_inf","HPD_sup")
    rownames(table) = names_x
    names(coeff) = names_x
  }
  list = list(table = table,betas = coeff)
  return(list)
}
