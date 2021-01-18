#' Check prior distributions
#'
#' @details The function is used to specify the prior-related arguments of the state-space modeling function \code{\link{run_ssm}}. Priors are specified in terms of distributions and associated parameters as implemented in the \pkg{rstan} package.
#' The available options are as follows:
#' \itemize{
#' \item \code{lognormal(mu,sigma)} (code = 1, required parameters = 2)
#' \item \code{normal(mu,sigma)} (code = 2, required parameters = 2)
#' \item \code{normal(mu,sigma)T(min,max)} (code = 201, required parameters = 4)
#' \item \code{chi_square(df)} (code = 3, required parameters = 1)
#' \item \code{inv_chi_square(df)} (code = 4, required parameters = 1)
#' \item \code{gamma(alpha,beta)} (code = 5, required parameters = 2)
#' \item \code{pareto(min,alpha)} (code = 6, required parameters = 2)
#' \item \code{uniform(min,max)} (code = 7, required parameters = 2)
#' }
#' This is an internal function, generally not to be called by the user.
#' @param priors (list) a list of arguments specifying priors for each parameter involved in the model (see Details). If \code{priors=NULL} then predefined priors will be returned.
#' @return a matrix contatining priors (numeric codes, see Details) and their parameters
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Define priors for all the paramaters 
#' priors_list <- list("lognormal(1,1)","normal(2,3)T(0,10)","normal(3,1)")
#' priors_out <- check_prior(priors_list)
#' print(priors_out)
#' 
#' ## Define priors for some of the paramaters 
#' priors_list <- list(NULL,"pareto(1,1.2)",NULL)
#' priors_out <- check_prior(priors_list)
#' print(priors_out)
#' 
#' ## Use pre-defined vague priors for all the parameters
#' priors_list <- list(NULL,NULL,NULL)
#' priors_out <- check_prior(priors_list)
#' print(priors_out)
#' }

check_prior <- function(priors=NULL){
  if(!is.list(priors))
     warning("Priors should be written as a list (see ?check_prior)")
  
  # lognormal=1, normal=2, normal_truncated=201, chisquare=3, inv_chi_square=4, gamma=5, pareto=6, uniform=7
  prior_avail <- c("lognormal","normal","chi_square","inv_chi_square","gamma","pareto","uniform")
  
  # Check NULL elements in the list and replace them with pre-defined priors
  iid <- which(sapply(priors,is.null)) 
  priors[iid] <- ifelse(iid==1,"lognormal(1,0.5)","normal(0,2)") #lognormal() is for the dummy intercept of the stimuli equation
  
  prior_dist <-mapply(function(i){
    strs <- strsplit(priors[[i]],"[,)T(]")[[1]][1]
  },1:length(priors),SIMPLIFY = FALSE)
  
  prior_params <- mapply(function(i){
    strs <- strsplit(priors[[i]],"[,)T(]")[[1]]
    params <- as.numeric(strs[2:length(strs)]); params <- params[!is.na(params)]
  },1:length(priors),SIMPLIFY = FALSE)
  
  priors_code <- match(prior_dist,prior_avail)
  
  # Check for normal_truncated by length of params
  iid <- which(sapply(prior_params,length)>2)
  priors_code[iid] <- 201
  
  prior_params_matrix <- matrix(0,length(prior_params),4)
  for(i in 1:length(prior_params)){
    x <- unlist(prior_params[[i]])
    prior_params_matrix[i,1:length(x)] <- x
  }
  
  prior_matrix <- cbind(priors_code,prior_params_matrix)
  
  return(prior_matrix)

}