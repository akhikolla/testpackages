

likeli_i <- function(yi, delta, X, Y){
  I <- nrow(delta)
  k <- ncol(delta)

    pi_i <- c()
    for(o in 1:I){
      pi_i[o] <- prob_pi(1,o,X,delta,Y,k)[yi[o]]
    }
    loglik <- prod(pi_i)
  
  
  loglik
  
}



estim_i<- function(x, yi, delta, Sigma){

    dens.norm <- dmvnorm(x, sigma = Sigma)
    likmat <- likeli_i(yi,delta,x[1],x[2])
    ret <- dens.norm * likmat * x

  return(ret)
}


factor_i <- function(x, yi, delta, Sigma){

  dens.norm <- dmvnorm(x, sigma = Sigma)
  likmat <- likeli_i(yi,delta,x[1],x[2])
  ret <- dens.norm * likmat 

  return(ret)
}


person_i <- function(person, model, Y, tol, limits, maxEval){

  delta <- model$delta
  Sigma <- model$Sigma
  
  yi <- Y[person,]
  
  est_i <- adaptIntegrate(estim_i,lowerLimit=-limits, upperLimit = limits, 
                             yi=yi, delta = delta, Sigma = Sigma, fDim = 2, tol = tol,
                          maxEval = maxEval)

est_i <- est_i$integral/adaptIntegrate(factor_i,lowerLimit=-limits, upperLimit = limits, 
                            yi=yi, delta = delta, Sigma = Sigma, tol = tol,
                            maxEval = maxEval)$integral
  
  est_i
}

#' Calculate Posterior Estimates for Person Parameters
#' 
#' Calculates posterior estimates for both person parameters, namely the ability parameters theta and the 
#' response style parameters gamma.
#' 
#' 
#' @param model Object of class \code{PCMRS}. 
#' @param cores Number of cores to be used in parallelized computation.
#' @param tol The maximum tolerance for numerical integration, default 1e-4. 
#' For more details see \code{\link{adaptIntegrate}}.
#' @param maxEval The maximum number of function evaluations needed in numerical integration.
#' If specified as 0 implies no limit. For more details see \code{\link{adaptIntegrate}}.
#' @param which Optional vector to specify that only for a subset of all persons the posterior estimate is calculated.
#' @return 
#' \item{}{Matrix containing all estimates of person parameters, both theta and gamma.} 
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://www.sg.tum.de/epidemiologie/team/schauberger/}
#' @seealso \code{\link{PCMRS}} \code{\link{PCMRS-package}}
#' @references Tutz, Gerhard, Schauberger, Gunther and Berger, Moritz (2018): 
#' Response Styles in the Partial Credit Model, \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/10.1177/0146621617748322}
#' @keywords PCMRS
#' @examples
#' \dontshow{
#' k <- 4; n <- 40; I <- 4
#' set.seed(1860)
#' Y <- as.data.frame(matrix(sample(1:k, I*n, TRUE),nrow = n))
#' Y <- data.frame(lapply(Y, as.ordered))
#' 
#' mini.ex <- PCMRS(Y, cores = 2)
#' 
#' persons <- person.posterior(mini.ex, cores = 1, which = 1)
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
person.posterior <- function(model, cores = 30, tol = 1e-4, maxEval = 600, which = NULL){

  Y <- model$Y

  if(!is.null(which)){
    Y <- Y[which,,drop = FALSE]
  }
  
  y_patterns <- apply(Y,1,paste,collapse="")
  y_pat_fac <- as.numeric(as.factor(y_patterns))

  
  Y_big <- matrix(as.numeric(as.matrix(Y)),ncol=ncol(Y))
  
  Y <- Y_big[!duplicated(y_patterns),,drop = FALSE]
  n <- nrow(Y)
  
  limits <- c(qnorm(0.9995,sd=sqrt(model$Sigma[1,1])),qnorm(0.9995,sd=sqrt(model$Sigma[2,2])))

  if(cores>1){
    cl <- makeCluster(cores,outfile="")
    
    clusterExport(cl, varlist = c("Y","model","factor_i",
                                  "estim_i","estim_i", "likeli_i","prob_pi",
                                  "prob_pir","lp","tol","limits","maxEval") 
                  ,envir = sys.frame(sys.nframe()))
    
    estimates <- parSapply(cl, seq(n), person_i, model = model, Y=Y, tol= tol,
                           limits = limits, maxEval = maxEval)
    stopCluster(cl)
  }else{
    estimates <- sapply(seq(n), person_i, model = model, Y=Y, tol= tol,
                        limits = limits, maxEval = maxEval)
  }
  
  new_order <- y_pat_fac[!duplicated(y_patterns)]
  estimates <- estimates[,order(new_order),drop = FALSE]
  estimates <- estimates[,y_pat_fac,drop = FALSE]

  rownames(estimates) <- c("theta", "gamma")
  colnames(estimates) <- rownames(model$Y)[which]
  
  t(estimates)
}