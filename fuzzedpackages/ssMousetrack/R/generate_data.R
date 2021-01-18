#' Generate datasets according to the model structure
#'
#' @details The function generates simulated datasets via Stan according to the model structure.
#' @param M (integer) number of simulated datasets
#' @param N (integer) length of the Y-trajectories
#' @param I (integer) number of individuals
#' @param J (integer) number of trials
#' @param K (array of integers) list of length Q of the number of levels for each categorical variable
#' @param Z.type (array of characters) list of length Q of the methods (symmetric or random) to generate the matrix (see \code{\link{generate_Z}})
#' @param Z.contrast (character) type of contrasts (default: treatment) for the model matrix Z (see \code{\link{model.matrix}})
#' @param Z.formula (character) a formula of the contrasts for the model matrix Z (see \code{\link{model.matrix}})
#' @param sigmax (numeric) fixed value for the model parameter sigmax
#' @param lambda (numeric) fixed value for the model parameter lambda
#' @param yT (numeric) position in angles of the target
#' @param yD (numeric) position in angles of the distractor
#' @param kappa_bnds (array) array containing the lower and upper bounds for the kappa parameter (\code{default = c(120,300)})
#' @param priors (list) a list of arguments specifying priors for each parameter involved in the model (see \code{\link{check_prior}}). If \code{priors="default"} then pre-defined priors will be used.
#' @param gfunction (character) type of link function between latent states and observed data: 'logistic', 'gompertz' (\code{default = 'logistic'}).  
#' @param ... other stan arguments (e.g., 'init', 'algorithm', 'sample_file'. See \code{\link[rstan:sampling]{sampling}}) 
#' @return a datalist containing simulated data and parameters
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Generate mouse-tracking data for an univariate experimental design 
#' ## with K = 3 categorical levels, J = 30 trials, I = 8 subjects
#' X1 <- generate_data(I=5,J=12,K=3,Z.formula="~Z1",M=50)
#' 
#' ## Generate mouse-tracking data for an univariate experimental design 
#' ## by varying priors of parameters
#' priors_list = list("normal(0,1)T(0,Inf)","normal(0,1)","normal(-2,0.5)")
#' X1 <- generate_data(I=5,J=12,K=3,Z.formula="~Z1",M=50,priors=priors_list)
#'
#' ## Generate mouse-tracking data with two experimental factors Z1 and Z2, J = 9 trials, 
#' ## K_Z1 = 3, K_Z2 = 3, I = 5 subjects
#' X2 <- generate_data(I=5,J=9,K=c(3,3),Z.formula="~Z1*Z2",
#' Z.type=c("symmetric","random"),M=50) # design with interaction
#' }

generate_data <- function(M=100,N=61,I=10,J=12,K=c(4),Z.type=c("symmetric"),Z.contrast="treatment",
                          Z.formula=NULL,sigmax=1,lambda=1,yT=pi/4,yD=(3*pi)/4,kappa_bnds=c(120,300),priors="default",gfunction=c("logistic","gompertz"),...){
  if(I<1| N<1 | J<1 | M<1)
    stop("Positive integers should be provided for I, J, N, M")
  if(length(K)<1)
    stop("At least an integer should be provided for K")
  if(length(Z.type)<1)
    stop("At least a method should be provided to generate Z")
  if(sigmax<0)
    stop("The parameter sigma must be greater then zero")
  if(lambda<0)
    stop("The parameter lambda must be greater then zero")
  if(yT<=0 | yT > pi)
    stop("The parameter yT must lie in the range (0,pi]")
  if(yD<=0 | yD > pi)
    stop("The parameter yD must lie in the range (0,pi]")
  if(yT >= yD)
    stop("The parameter yD must be greater yT")
  if(is.null(Z.formula))
    stop("The model.matrix formula for Z must be provided")
  
  gfunction=match.arg(gfunction)

  X.design <- generate_design(I,J,K,Z.type)
  Z.matrix <- stats::model.matrix(formula(Z.formula),data = X.design,contrasts.arg = Z.contrast)
  
  if(length(priors)==1 && priors=="default"){
    priors <- rep(list(NULL),dim(Z.matrix)[2])}

  lb <- 0.1 #to avoid the case y[n]~0

  datastan <- list(
    I = I,
    N = N,
    J = J,
    KK = dim(Z.matrix)[2],
    Y = array(0,c(N,I*J)),
    sigmaz = rep(sigmax,I),
    bnds = matrix(1,I*J,1)%*%matrix(c(lb,pi-lb,(pi-lb)-lb),1,3),
    D = Z.matrix,
    lambda_vec = rep(lambda,I*J),
    priors_matrix = check_prior(priors),
    pT = yT,
    pD = yD,
    pC = (yT+yD)/2,
    kappa_lb = kappa_bnds[1],
    kappa_ub = kappa_bnds[2]
  )
  
  if(gfunction=="logistic"){
    out <- rstan::sampling(stanmodels$simulate_data_log, data = datastan,chain=1,iter=M,warmup=M/4,...)
  }else if(gfunction=="gompertz"){
    out <- rstan::sampling(stanmodels$simulate_data_gomp, data = datastan,chain=1,iter=M,warmup=M/4,...)
  }
  
  data_out <- rstan::extract(out,pars=c("b","z","mu","y_sim","gamma","dy_sim"))
  
  datasim <- list(
    I = I,
    N = N,
    J = J,
    K = K,
    Gfunction = gfunction,
    params = list(sigmax=rep(sigmax,I),lambda=rep(lambda,I*J),gamma=data_out$gamma,beta=data_out$b),
    data = list(Y=data_out$y_sim,X=data_out$z,MU=data_out$mu,D=data_out$dy_sim,Z=Z.matrix),
    design = X.design
  )

  return(datasim)
}
