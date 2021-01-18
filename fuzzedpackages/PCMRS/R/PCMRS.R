


## response function for acat
responseFun <- function(eta){
  q <- length(eta)
  eta.help <- matrix(rep(c(0,eta),each=q+1),ncol=q+1)
  eta.help[upper.tri(eta.help)] <- 0
  pi <- cumprod(c(1,exp(eta[-q])))/sum(apply(exp(eta.help),1,prod))
  pi
}

## create responses for acat from ordinal values
createResponse <- function(Y){
  c(t(model.matrix(~0+Y)[,-length(levels(Y))]))
}

#' Model Response Styles in Partial Credit Models
#' 
#' Performs PCMRS, a method to model response styles in Partial Credit Models
#' 
#' 
#' @param Y Data frame containing the ordinal item response data (as ordered factors), one row per obeservation, one column per item.
#' @param Q Number of nodes to be used (per dimension) in two-dimensional Gauss-Hermite-Quadrature.
#' @param scaled Should the scaled version of the response style parameterization be used? Default is \code{TRUE}.
#' @param method Specifies optimization algorithm used by \code{\link{optim}}, either 
#' \code{L-BFGS-B} or \code{nlminb}.
#' @param cores Number of cores to be used in parallelized computation.
#' @return 
#' \item{delta}{Matrix containing all item parameters for the PCMRS model, one row
#' per item, one column per category.} 
#' \item{Sigma}{2*2 covariance matrix for both random effects, namely the ability parameters theta and the
#' response style parameters gamma.}
#' \item{delta.PCM}{Matrix containing all item parameters for the simple PCM model, one row
#' per item, one column per category.} 
#' \item{sigma.PCM}{Estimate for variance of ability parameters theta in the simple PCM model.}
#' \item{Y}{Data frame containing the ordinal item response data, one row per obeservation, one column per item.} 
#' \item{scaled}{Logical, \code{TRUE} if scaled version of the response style parameterization is used.} 
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://www.sg.tum.de/epidemiologie/team/schauberger/}
#' @seealso \code{\link{person.posterior}} \code{\link{PCMRS-package}}
#' @references Tutz, Gerhard, Schauberger, Gunther and Berger, Moritz (2018): 
#' Response Styles in the Partial Credit Model, \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/10.1177/0146621617748322}
#' @keywords PCMRS
#' @examples
#' \dontshow{
#' k <- 4; n <- 50; I <- 4
#' set.seed(1860)
#' Y <- as.data.frame(matrix(sample(1:k, I*n, TRUE),nrow = n))
#' Y <- data.frame(lapply(Y, as.ordered))
#' 
#' mini.ex <- PCMRS(Y, cores = 2)
#' mini.ex
#' }
#' \dontrun{
#' ################################################
#' ## Small example to illustrate model and person estimation
#' ################################################
#' 
#' data(tenseness)
#' 
#' set.seed(5)
#' samples <- sample(1:nrow(tenseness), 100)
#' tense_small <- tenseness[samples,1:4]
#' 

#' m_small <- PCMRS(tense_small, cores = 2)
#' m_small
#' plot(m_small)
#' 
#' persons <- person.posterior(m_small, cores = 2)
#' plot(jitter(persons, 100))
#' 
#' ################################################
#' ## Example from Tutz et al. 2017:
#' ################################################
#' 
#' data(emotion)
#' m.emotion <- PCMRS(emotion)
#' m.emotion
#' 
#' plot(m.emotion)
#' }
PCMRS <- function(Y, Q = 10, scaled = TRUE, method = c("L-BFGS-B", "nlminb"), cores = 30){

  ## For cran version X is eliminated from command and set NULL automatically!
  X <- NULL
  
  if(!is.data.frame(Y))
    stop("Y has to be a data.frame!")
  
  all.factors <- TRUE
  for(ii in 1:ncol(Y)){
    all.factors <- all(all.factors, is.ordered(Y[,ii]))
  }
  if(!all.factors)
    stop("All variables in Y have to be ordered factors!")
  
  
  method <- match.arg(method)

  ## initalize response vector and other parameters
  y.vec <- as.factor(c(t(Y)))
  k <- length(levels(y.vec))
  q <- k-1
  n <- nrow(Y)
  I <- ncol(Y)
  
  if(is.null(X)){
    X <- matrix(0,0,0)
  }
  pX <- ncol(X)
  
  ## create response for acat family
  response <- createResponse(y.vec)

  ## initialize starting values, 
  ## just use a partial credit model for delta
  ## and the sigma for both sigmas,
  m.ltm <- gpcm(Y, constraint = "1PL")
  coef.ltm <- coef(m.ltm)
  sigma1.start <- coef.ltm[1,q+1]
  delta.start <- c(t(coef.ltm[,-(q+1)]))*sigma1.start
  beta.start <- rep(0,2*pX)
  alpha.start <- c(delta.start, sigma1.start^2, 0, 0.5, beta.start)
  
  ## get nodes and weights
  her_poly <- gauss.quad(Q, "hermite")
  nodes <- her_poly$nodes
  weights <- her_poly$weights * exp(nodes^2)
  weights <- weights %*% t(weights)
  node.probs = dnorm(nodes)%*%t(dnorm(nodes))
  

## lower bounds for variance parameters
## upper and lower bound for correlation
  l.bound <- rep(-Inf,length(alpha.start))
  l.bound[c(length(delta.start)+1,length(delta.start)+3)] <- 1e-10
  l.bound[length(delta.start)+2] <- -1+1e-4
  u.bound <- rep(Inf,length(alpha.start))
  u.bound[length(delta.start)+2] <- 1-1e-4
  

## change scaling for variance and correlation parameters
 par.scale <- rep(1,length(alpha.start))

########################
if(method == "L-BFGS-B"){
  m.opt <- try(optim(par = alpha.start, fn = loglikEPCM, gr = scoreEPCM,
                  Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
                  pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes,
                  scaled = as.numeric(scaled),
                  cores = cores, lower = l.bound, upper = u.bound, method="L-BFGS-B",
                  control = list(parscale = par.scale)))

  if( (class(m.opt) == "try-error")){
    par.scale <- rep(0.1,length(alpha.start))
    
    m.opt <- try(optim(par = alpha.start, fn = loglikEPCM, gr = scoreEPCM, 
                       Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
                       pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes, 
                       scaled = as.numeric(scaled),
                       cores = cores, lower = l.bound, upper = u.bound, method="L-BFGS-B",
                       control = list(parscale = par.scale)))
 }
 
 if( (class(m.opt) == "try-error")){
   par.scale <- rep(1,length(alpha.start))
   method <- "nlminb"
 }
}

if(method == "nlminb"){
  m.opt <- try(nlminb(start = alpha.start, objective = loglikEPCM, gradient = scoreEPCM, 
                  Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
                  pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes, 
                  scaled = as.numeric(scaled),
                  cores = cores, lower = l.bound, upper = u.bound, scale= par.scale))
  
  if( (class(m.opt) == "try-error")){
    par.scale <- rep(0.1,length(alpha.start))
    
    m.opt <- try(nlminb(start = alpha.start, objective = loglikEPCM, gradient = scoreEPCM, 
                    Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
                    pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes, 
                    scaled = as.numeric(scaled),
                    cores = cores, lower = l.bound, upper = u.bound, scale= par.scale))

  }
  if( (class(m.opt) == "try-error")){
    par.scale <- rep(0.01,length(alpha.start))
    
    m.opt <- try(nlminb(start = alpha.start, objective = loglikEPCM, gradient = scoreEPCM, 
                        Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
                        pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes, 
                        scaled = as.numeric(scaled),
                        cores = cores, lower = l.bound, upper = u.bound, scale= par.scale))
    
  }
}

########################
  
delta <- Sigma <- NA

if( (class(m.opt) != "try-error")){
  ## extract results and prepare return
  coefs <- m.opt$par
# browser()
#   n0 <- grad(loglikEPCM, coefs,
#              Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
#              pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes,
#              scaled = as.numeric(scaled),
#              cores = 1)
# 
#   a0 <- scoreEPCM(alpha =coefs,
#                   Q = Q, q = q, I = I, n = n, Y = response, X = X, pall = length(alpha.start),
#                   pX = pX, GHprobs = node.probs, GHweights = weights, GHnodes = nodes,
#                   scaled = as.numeric(scaled),
#                   cores = 1)
# 
#   print( rbind(c(a0),c(n0)))
  
  delta <- matrix(coefs[1:(length(coefs)-3-2*pX)], byrow=TRUE, ncol = q)
  colnames(delta) <- paste("Catgr",1:q,sep=".")
  rownames(delta) <- colnames(Y)
  
  sigma <- coefs[(length(coefs)-2-2*pX):(length(coefs)-2*pX)]
  if(sigma[3]==0){
    Sigma <- matrix(c(sigma[1],0,0,0), ncol = 2)
  }else{
    Sigma <- matrix(c(sigma[1], sigma[2]*sqrt(sigma[1])*sqrt(sigma[3]),
                      sigma[2]*sqrt(sigma[1])*sqrt(sigma[3]), sigma[3]),ncol=2)
  }
}

colnames(Sigma) <- rownames(Sigma) <- c("theta","gamma")

  delta.PCM <- matrix(delta.start, ncol=q, byrow=TRUE)
  colnames(delta.PCM) <- paste("Catgr",1:q,sep=".")
  rownames(delta.PCM) <- colnames(Y)

  sigma.PCM <- sigma1.start^2
  
  ## define parameters of covariate influence on trait parameters
  beta.theta <- beta.gamma <- NULL
  if(pX>0){
    beta.theta <- coefs[(length(coefs)-2*pX+1):(length(coefs)-pX)]
    beta.gamma <- coefs[(length(coefs)-pX+1):(length(coefs))]
  }

  ret.list <- list(delta = delta, Sigma = Sigma, beta.theta = beta.theta, 
                   beta.gamma = beta.gamma, delta.PCM = delta.PCM, sigma.PCM = sigma.PCM,
                   Y = Y, scaled = scaled)

  class(ret.list) <- "PCMRS"
  
  return(ret.list)
}
