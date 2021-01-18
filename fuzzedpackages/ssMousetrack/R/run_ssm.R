#' State-space modeling of mouse-tracking trajectories via Stan 
#'
#' @details The function drawns samples from the posterior distribution of the model parameters. Note that, the current version of \pkg{ssMousetrack} package requires the number of stimuli J to be the same over the subjects \eqn{i=1,...,I}.
#' @param N (integer) length of the Y-trajectories
#' @param I (integer) number of individuals
#' @param J (integer) number of trials
#' @param Y (matrix) N x JI matrix of observed trajectories
#' @param D (matrix) N x JI matrix of delta values for the observed trajectories
#' @param Z (matrix) matrix of contrasts associated to the experimental design (see \code{\link{generate_design}})
#' @param sigmax (numeric) fixed value for the model parameter sigmax
#' @param lambda (numeric) fixed value for the model parameter lambda
#' @param y_T (numeric) position in angles of the target
#' @param y_D (numeric) position in angles of the distractor
#' @param priors (list) a list of arguments specifying priors for each parameter involved in the model (see \code{\link{check_prior}}). If \code{priors="default"} then pre-defined tpriors will be used.
#' @param gfunction (character) type of link function between latent states and observed data: 'logistic', 'gompertz' (\code{default = 'logistic'}).  
#' @param kappa_bnds (array) array containing the lower and upper bounds for the kappa parameter (\code{default = c(5,300)})
#' @param nchains (integer) number of chains for the MCMC algorithm
#' @param niter (integer) number of iterations for each chain
#' @param nwarmup (integer) number of warmup/burnin iterations per chain
#' @param ncores (integer) number of cores to use when executing the chains in parallel. The default option ncores="AUTO" will automatically determine the number of cores via the \code{parallel} package
#' @param stan_object (boolean) if \code{stan_object=TRUE}, the object of S4 class stanfit representing the fitted results will be saved as stan_object.rda
#' @param ... other stan arguments (e.g., 'init', 'algorithm', 'sample_file'. See \code{\link[rstan:sampling]{sampling}})
#' @return a datalist containing the posterior samples for the model parameters along with the main Stan output 
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Fit a state-space model using simulated data 
#' # Generate mouse-tracking data for an univariate experimental design with K = 3 categorical levels, 
#' # J = 12 trials, I = 5 subjects
#' X1 <- generate_data(I=5,J=12,K=3,Z.formula="~Z1")
#' iid <- 23 # keep just one dataset from the simulated set of datasets
#' # Run the state-space model on the chosen dataset
#' X1_fit <- run_ssm(N = X1$N,I = X1$I,J = X1$J,Y = X1$data$Y[iid,,],D = X1$data$D[iid,,],
#' Z = X1$data$Z)
#' 
#' ## Fit a state-space model using the experimental dataset language
#' # The dataset is ready to be used and it does not need to be pre-processed (preprocess=FALSE). 
#' # In this case, the function prepare_data just computes the observed radians from 
#' # the x-y trajectories
#' X2 <- prepare_data(X = language, preprocess = FALSE, Z.formula = "~condition")
#' # Run the state-space model on the chosen dataset
#' X2_fit <- run_ssm(N = X2$N,I = X2$I,J = X2$J,Y = X2$Y,D = X2$D,Z = X2$Z,
#' niter=5000,nchains=2)
#' 
#' ## Fit a state-space model using the experimental dataset congruency
#' # The dataset needs to be pre-processed (preprocess=TRUE)
#' X3 <- prepare_data(X = congruency, preprocess = TRUE, 
#' Z.formula = "~congruency+plausibility") # additive design
#' # Define priors of the model parameters 
#' KK <- dim(X3$Z)[2] # number of model parameters implied by the design matrix Z
#' priors_list <- list("lognormal(1,0.5)","pareto(3,5.25)","normal(0,2.5)")
#' # note that length(priors_list) = KK 
#' # Run the state-space model on the chosen dataset
#' X3_fit <- run_ssm(N = X3$N,I = X3$I,J = X3$J,Y = X3$Y,D = X3$D,Z = X3$Z,
#' niter=10000,nwarmup=3500,priors=priors_list,nchains=4)
#' }


run_ssm <- function(N,I,J,Y=NULL,D=NULL,Z=NULL,sigmax=1,lambda=1,y_T=pi/4,y_D=(3*pi)/4,priors="default",gfunction=c("logistic","gompertz"),kappa_bnds=c(5,300),nchains=1,niter=2000,nwarmup=500,ncores="AUTO",stan_object=FALSE,...){
  if(I<1| N<1 | J<1)
    stop("Positive integers should be provided for I, J, N, M")
  if(is.null(Z))
    stop("The design matrix Z must be provided")
  if(is.null(Y))
    stop("The observed matrix Y must be provided")
  if(is.null(D))
    stop("The matrix of distances D must be provided")
  if(ncores=="AUTO")
    ncores <- parallel::detectCores()
  if(length(priors)==1 && priors=="default")
    priors <- rep(list(NULL),dim(Z)[2])
  
  gfunction=match.arg(gfunction)
  
  lb <- 0.1 #to avoid the case y[n]~0
  
  datastan <- list(
    I = I,
    N = N,
    J = J,
    KK = dim(Z)[2],
    Y = Y,
    DY = D,
    sigmaz = rep(sigmax,I),
    bnds = matrix(1,I*J,1)%*%matrix(c(lb,pi-lb,(pi-lb)-lb),1,3),
    D = Z,
    a = rep(1,I*J),
    lambda_vec = rep(lambda,I*J),
    Am = kronecker(diag(1,I),rep(1/J,J)),
    kappa_lb = kappa_bnds[1],
    kappa_ub = kappa_bnds[2],
    priors_matrix = check_prior(priors)
  )

  if(gfunction=="logistic"){
    out <- rstan::sampling(stanmodels$fit_model_log,chains=nchains,iter=niter,warmup=nwarmup,data=datastan,cores=ncores,...)
  }else if(gfunction=="gompertz"){
    out <- rstan::sampling(stanmodels$fit_model_gomp,chains=nchains,iter=niter,warmup=nwarmup,data=datastan,cores=ncores,...)
  }
    
  data_out <- rstan::extract(out,pars=c("b","z_pred","y_star","gamma","z_s_upd"))

  gamma_out = data.frame(data_out$gamma)
  names(gamma_out) = paste("gamma",seq(1,dim(gamma_out)[2]),sep="")
  
  datafit <- list(
    I = I,
    N = N,
    J = J,
    Gfunction = gfunction,
    params = list(sigmax=sigmax,lambda=lambda,kappa_bnds=kappa_bnds,gamma=gamma_out,beta=data_out$b,kappa),
    data = list(Y=Y,X=data_out$z_pred,MU=data_out$y_star,D=D,Z=Z,X_smooth=data_out$z_s_upd),
    stan_table = rstan::monitor(out)
  )
  
  if(stan_object==TRUE){
    save(out,file = paste(getwd(),"/stan_object.rda",sep=""),compress = TRUE)
    }
  
  return(datafit)
}
