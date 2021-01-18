

# draw one sample of chpts and states from the posterior
draw.from.post <- function(logprobs, cpt.locs, types, params){

  n <- length(logprobs)
  l.probs <- logprobs[[n]]
  seg.locs <- cpt.locs[[n]]
  seg.types <- types[[n]]

  ind <- 1:length(l.probs)
  a <- sample( ind , size=1 , replace = TRUE, exp( l.probs ) )

  # state amd location
  curr.state <- seg.types[a]
  t <- seg.locs[a]

  STATES <- curr.state
  DRAW <- t

  k_N = params[1]
  p_N = params[2]
  k_A = params[3]
  p_A = params[4]
  pi_N = params[5]
  pi_A = 1 - pi_N

  while (t > 1){

    if(curr.state==0){

      l.probs <- logprobs[[t]][ types[[t]] == 1 ]
      i <- cpt.locs[[t]][ types[[t]] ==1 ]
      back.dens <- log(pi_N) + dnbinom(t-i-1,k_A,p_A,log=TRUE) - pnbinom(t-i-2,k_A,p_A, lower.tail=FALSE , log.p=TRUE)
      l.probs <- l.probs + back.dens

      ind <- 1:length(l.probs)
      a <- sample( ind , size=1 , replace = TRUE, exp( l.probs ) )

      # state amd location
      curr.state <- 1
      t <- i[a]
      STATES <- c(STATES,curr.state)
      DRAW <- c(DRAW,t)
    }
    else{

      # currently in abnormal state poss transitions to normal or abnormal
      i.N <- cpt.locs[[t]][ types[[t]] == 0 ]
      i.A <- cpt.locs[[t]][ types[[t]] == 1 ]

      l.probsA <- logprobs[[t]][  types[[t]] == 1 ]
      back.densA <- log(pi_A) + dnbinom(t-i.A-1,k_A,p_A,log=TRUE) - pnbinom(t-i.A-2,k_A,p_A, lower.tail=FALSE , log.p=TRUE)
      l.probsA <- l.probsA + back.densA
      # normal points
      l.probsN <- logprobs[[t]][  types[[t]] == 0 ]
      back.densN <- dnbinom( t-i.N-1 , k_N , p_N , log=TRUE ) - pnbinom( t-i.N-2 , k_N , p_N , lower.tail=FALSE , log.p=TRUE )
      l.probsN <- l.probsN + back.densN
      # now normalise
      i.all <- c( i.N , i.A )
      t.all <- c( rep(0,length(i.N)) , rep(1,length(i.A)) )
      l.probs <- c( l.probsN , l.probsA )
      c <- max(l.probs)
      l.probs.norm <- l.probs - ( c + log( sum( exp( l.probs - c ) ) ) )

      ind <- 1:length(l.probs.norm)
      a <- sample( ind , size=1 , replace = TRUE, exp( l.probs.norm ) )

      # state amd location
      curr.state <- t.all[a]
      t <- i.all[a]
      STATES <- c(STATES,curr.state)
      DRAW <- c(DRAW,t)
    }

  }

  DRAW[DRAW==0] <- 1
  newlist <- list("draws"=c(n,DRAW),"states"=STATES)
  return(newlist)

}

# loss function
loss <- function(gamma, probvec){

  p.gamma <- 1/(1+gamma)
  probvec[probvec < p.gamma ] <- 0
  probvec[probvec>0] <- 1
  return(probvec)

}


get.states <- function(segs, states, n){

  state.vec <- numeric(n)
  state.vec[ segs[1]:segs[2] ] <- states[1]
  for (i in 2:( length(segs) - 1 ) ){
    state.vec[ ( segs[i]-1 ):segs[i+1] ] <- states[i]
  }

  return(state.vec)

}


Rsampler <- function(res, gamma = 1/3, no.draws=1000){
  
  logprobs <- res$weights
  cpt.locs <- res$locations
  types <- res$type
  params <- res$params
  
  # sampler params
  mu_seq = res$sampler_params$mu_seq
  N = res$sampler_params$N
  S = res$sampler_params$S
  p = res$sampler_params$p
  
  n = length(logprobs)
  prob.states <- matrix(nrow=no.draws, ncol=n)
  # for the heatmap  
  sampled.res = list()
  for (i in 1:no.draws){
    
    f <- draw.from.post(logprobs, cpt.locs, types, params)
    segs <- f$draws
    states <- f$states
    
    # for marginal probabilities
    prob.states[i,] <- get.states(segs, states, n)
    
    # for the heatmap
    end <- segs[which(states == 1)]
    begin <- segs[which(states == 1) + 1]
    tmpmat <- matrix(nrow = 2, ncol = length(end))
    tmpmat[1,] <- begin
    tmpmat[2,] <- end
    sampled.res[[i]] <- tmpmat
    
    
  }
  
  margprob = apply(prob.states,2,sum)/no.draws
  segmentation = loss(gamma, margprob)
  
  welem = which(segmentation == 1)
  if (length(welem) > 0){
    
    # put into data frame
    df = data.frame("start"=NA, "end"=NA, "LogMargLike"=NA)
    wdiff = c(0, which(diff(welem) != 1), length(welem))
    for (i in 1:(length(wdiff)-1) ){
      start = welem[(wdiff[i] + 1)]
      end = welem[wdiff[(i+1)]]
      ml = P.A(start, end, mu_seq, N, S, p)
      tempdf = data.frame("start"=start,
                          "end"=end,
                          "LogMargLike"=ml)
      df = rbind(df, tempdf)
    }
    df = df[-1,]
    df = df[order(df$LogMargLike, decreasing = T), ]
    return(list(df,margprob,sampled.res))
    
  }
  else{
    df = data.frame()
    return(list(df,margprob,sampled.res))
  }
  
}





# Resampling function -- takes a vector of weights that are < alpha 
# to resample and does SRC resampling, then returns a vector with resampled 
# weights 

resample <- function(to.resample,alpha){ 
  
  log.alpha <- log(alpha)
  log.u <- log.alpha + log( runif(1) )
  k <- 1

  while (k <= length(to.resample) ){
    
    # want to find u <- u - w then if u <= 0 
    # u <- u + alpha
    # IN OUR CASE IF u <= w
    # find u <- u + alpha - w
    # as working with logs
    # SO first check whether log(u) < log(w) 
  
    if ( log.u < to.resample[k] ){
      # u <- u + alpha - w
      # first find log(alpha-w) label as log.alpha.weight
      temp = c( log.alpha ,  to.resample[k] )
      c = max(temp)
      log.alpha.weight <- c + log( exp( temp[1] - c ) - exp( temp[2] - c )  )
      temp = c( log.u , log.alpha.weight )
      c = max(temp)
      log.u <- c + log( sum( exp( temp - c ) ) )
      to.resample[k] <- log.alpha
    }
  
    else{
      # u <- u - w
      # and w is not resampled, given weight 0 (-Inf for log(w))
      temp <- c( log.u , to.resample[k] )
      c = max(temp)
      log.u <- c + log( exp( temp[1] - c ) - exp( temp[2] - c )  )
      to.resample[k] <- - Inf
    
    }
  
    k <- k+1
    
  }

  return(to.resample)
  
}


# Calculating log P.A(s,t) for segment (s,t) integrated on uniform prior over mu_seq
# done without normal part at front
P.A <- function(s, t, mu_seq, N, S, p){

  # prior value mu_dens
  mu_dens <- 1/( tail(mu_seq,1) - mu_seq[1] )

  # width of rectangle
  mu_wid <- diff(mu_seq)[1]

  vec <- numeric(length(mu_seq))

  # evaluate at each point of grid
  # do this as typically smaller than dimension
  # evaluating log of quantity
  for (k in 1:length(mu_seq)){

    Z = mu_seq[k] * (S[(t+1),] - S[s,] - mu_seq[k] * (t-s+1)/2) + log(p) - log(1-p)
    Q = log( 1 + exp(Z) )
    wQ = which(Q == Inf)
    Q[wQ] = Z[wQ]
    vec[k] <- N*log(1-p) + sum(Q)
    
  }

  # finding sum of logs -- for numerical instability
  cmax <- max( vec )

  marg.like <- cmax + log( sum( exp( vec - cmax ) ) ) + log(mu_dens) + log(mu_wid)

  return(marg.like)

}


Rbard <- function(data, bardparams, mu_seq, alpha = 1e-4){

  if (length(bardparams) != 6){
    stop("Not enough params should be a vector of length 6.")
  }

  k_N = bardparams[1]
  p_N = bardparams[2]
  k_A = bardparams[3]
  p_A = bardparams[4]
  pi_N = bardparams[5]
  affected_dim = bardparams[6]
  pi_A = 1 - pi_N

  ## stat distribution qN , qA
  EN <- ( k_N * (1-p_N) )/p_N
  EA <- ( k_A * (1-p_A) )/p_A

  ldenom <- log( pi_N*EN + EA )
  qA <- log(pi_N) + log(EN) - ldenom
  qN <- log(EA) - ldenom

  # length and dimension of data
  n = dim(data)[1]
  N = dim(data)[2]
  # data summaries etc
  S <- rbind( rep(0,dim(data)[2]) , apply(data , 2 , cumsum) )
  S_2 <- cumsum( c( 0 , rowSums(data^2) ) )
  # p is used to calc the marginal like for abnormal
  # ratio of abnormal profiles
  p <- affected_dim/N

  # useful to define log of resampling probability
  log.alpha <- log(alpha)
  # lists for filtering, probs, location and types of segment
  weights <- vector("list",n)
  locations <- vector("list",n)
  type <- vector("list",n)
  # type - 0 for normal segment 1 for abnormal segment

  # initial N
  curr.locations <- c()
  curr.weights <- c()
  curr.type <- c()

  curr.locations[1] <- 0
  curr.weights[1] <- qN - 0.5 * ( S_2[2] - S_2[1] )
  curr.type[1] <- 0

  # initial A
  curr.locations[2] <- 0
  curr.weights[2] <- qA - 0.5 * ( S_2[2] - S_2[1] ) + P.A(1, 1, mu_seq, N, S, p)
  curr.type[2] <- 1

  c <- max(curr.weights)
  weights[[1]] <- curr.weights - ( c + log( sum( exp( curr.weights - c ) )  ) )
  locations[[1]] <- curr.locations
  type[[1]] <- curr.type

  for (t in 2:n){

    prev.type <- type[[t-1]]
    prev.locs <- locations[[t-1]]
    prev.weights <- weights[[t-1]]

    curr.type <- type[[t]]
    curr.locs <- locations[[t]]
    curr.weights <- weights[[t]]

    # label support points makes it easier to fnd positions in matrices
    sup.point.N <- which( prev.type == 0 )
    sup.point.A <- which( prev.type == 1 )

    # propagate normal particles
    # carry on in N state, nbinom.len
    # and P_N marginal likelihood

    i <- prev.locs[sup.point.N]
    nbinom.len <- pnbinom(t-i-1,k_N,p_N,log.p=T,lower.tail=F) - pnbinom(t-i-2,k_N,p_N,log.p=T,lower.tail=F)
    curr.weights[sup.point.N] <- prev.weights[sup.point.N] + nbinom.len - 0.5 * ( S_2[ t+1 ] - S_2[ t ] )

    # propagate abnormal particles
    # carry on in A state, nbinom.len
    # Marginal likelihood PA calculated at i+1:t and i+1:t-1
    for (j in sup.point.A){

      i <- prev.locs[j]
      nbinom.len <- pnbinom(t-i-1,k_A,p_A,log.p=T,lower.tail=F) - pnbinom(t-i-2,k_A,p_A,log.p=T,lower.tail=F)
      curr.weights[j] <- prev.weights[j] + nbinom.len - 0.5 * ( S_2[ t+1 ] - S_2[ t ] ) + P.A(i+1, t, mu_seq, N, S, p) - P.A(i+1, t-1, mu_seq, N, S, p)

    }

    ################################################
    # Calc other point chpt - C_t = t-1 , B_t = N/A
    ################################################

    ####
    # B_t = N
    ####
    # transition from abnormal to normal/ if there are no abnormal pts
    # give it a weight of zero i.e. log weight = - Inf
    i <- prev.locs[sup.point.A]

    if( length(i)==0 ){
      C_N <- - Inf
    }
    else{
      # weight to propogate
      W <- prev.weights[sup.point.A] + dnbinom(t-i-1,k_N,p_N,log=T) - pnbinom(t-i-2,k_N,p_N,log.p=T,lower.tail=F) + log(pi_N)

      # transition from A - N
      # (only trans possible to get to N)
      # prob of pi.N
      cmax <- max( W )
      C_N <- - 0.5*( S_2[ t+1 ] - S_2[ t  ] ) + cmax + log( sum( exp(W - cmax) ) )
    }

    ####
    # B_t = A
    ####
    # transition to abnormal, can be A-A or N-A
    i.N <- prev.locs[sup.point.N]
    i.A <- prev.locs[sup.point.A]

    if ( length(i.N) == 0 ){

      # no normal particles just use A particles
      W <- prev.weights[sup.point.A] + dnbinom( t - i.A -1  ,k_A,p_A,log=T) - pnbinom(t - i.A - 2,k_A,p_A,log.p=T,lower.tail=F)
      cmax <- max( W )
      C_A <- P.A(t, t, mu_seq, N, S, p) - 0.5*( S_2[ (t+1) ] - S_2[ t  ] ) + log(pi_A) + cmax + log( sum( exp(W - cmax) ) )

    }

    else if ( length(i.A) == 0 ) {

      # no abnormal particles just use N particles
      W <- prev.weights[sup.point.N] + dnbinom( t - i.N -1 ,k_N,p_N,log=T) - pnbinom(t - i.N - 2,k_N,p_N,log.p=T,lower.tail=F)
      cmax <- max( W )
      C_A <- P.A(t, t, mu_seq, N, S, p) - 0.5*( S_2[ (t+1) ] - S_2[ t  ] ) + cmax + log( sum( exp(W - cmax) ) )

    }

    else {

      # each have some N and A particles do seperately
      W <- prev.weights[sup.point.A] + dnbinom( t - i.A-1  ,k_A,p_A,log=T) - pnbinom(t - i.A - 2 ,k_A,p_A,log.p=T,lower.tail=F)
      cmax <- max( W )
      ABS.PART <- log(pi_A) + cmax + log( sum( exp(W - cmax) ) )
      W <- prev.weights[sup.point.N] + dnbinom( t - i.N -1  ,k_N,p_N,log=T) - pnbinom(t - i.N - 2 ,k_N,p_N,log.p=T,lower.tail=F)
      cmax <- max( W )
      NORM.PART <- cmax + log( sum( exp(W - cmax) ) )
      # put these together
      part <- c( NORM.PART , ABS.PART )
      cmax <- max(part)
      C_A <- P.A(t, t, mu_seq, N, S, p) - 0.5*( S_2[ (t+1) ] - S_2[ t  ] ) + cmax + log( sum( exp( part - cmax ) ) )

    }

    ###############################################
     ## resample if neccessary only when t > 20 ##
    ###############################################

    if ( t > 20 ){

      temp.weights <- c( curr.weights , C_N , C_A )
      temp.locs <- c( prev.locs , (t-1) , (t-1) )
      temp.type <- c( prev.type ,  0 , 1 )

      ## log normalized weights ##
      c <- max( temp.weights )
      lognorm.weights <-  temp.weights - ( c + log( sum(exp(temp.weights-c)) ) )

      ##############################
      # stratified resampling part #
      ##############################
      # take all the log weights that are < alpha and resample them
      to.resample <- lognorm.weights[ lognorm.weights < log.alpha ]
      # resampling function in resample.r
      lognorm.weights[ lognorm.weights < log.alpha ] <-  resample( to.resample , alpha )

      # remove weights for which = -Inf
      lognorm.weights[lognorm.weights == -Inf] <- NA
      updated.locs <- temp.locs[ !is.na(lognorm.weights) ]
      updated.type <- temp.type[ !is.na(lognorm.weights) ]

      # renormalise weights - is this necessary?
      temp.weights <- lognorm.weights[!is.na(lognorm.weights)]
      c <- max( temp.weights )
      renorm.weights <-  temp.weights - ( c + log( sum(exp(temp.weights-c)) ) )

      # put back in vectors (ordered)
      weights[[t]] <- renorm.weights
      locations[[t]] <- updated.locs
      type[[t]] <- updated.type

    }

    #########################################
         ## No resampling t <= 20 ##
    #########################################

    else{

      temp.weights <- c( curr.weights , C_N , C_A )

      # find log normalised weights #
      c <- max( temp.weights )
      lognorm.weights <-  temp.weights - ( c + log( sum(exp(temp.weights-c)) ) )

      weights[[t]] <- lognorm.weights
      locations[[t]] <- c( prev.locs , (t-1) , (t-1) )
      type[[t]] <- c( prev.type , 0 , 1  )

    }

  }

  sampler_params = list("mu_seq"=mu_seq, "N"=N, "S"=S, "p"=p)
  newList <- list("weights" = weights,
                  "locations" = locations,
                  "type" = type,
                  "params" = bardparams[1:5],
                  "sampler_params" = sampler_params)
  return(newList)

}



.bard.class<-setClass("bard.class",representation(data="array",
                                                  p_N="numeric",
                                                  p_A="numeric",
                                                  k_N="numeric",
                                                  k_A="numeric",
                                                  pi_N="numeric",
                                                  alpha="numeric",
                                                  paffected="numeric",
                                                  lower="numeric",
                                                  upper="numeric",
                                                  h="numeric",
                                                  Rs="list"))



bard.class<-function(data,p_N,p_A,k_N,k_A,pi_N,alpha,paffected,lower,upper,h,Rs,...)
{
    .bard.class(data=data,
                p_N=p_N,
                p_A=p_A,
                k_N=k_N,
                k_A=k_A,
                pi_N=pi_N,
                alpha=alpha,
                paffected=paffected,
                lower=lower,
                upper=upper,
                h=h,
                Rs=Rs)
}



.bard.sampler.class<-setClass("bard.sampler.class",
                              representation(bard.result="bard.class",
                              gamma="numeric",
                              num_draws="numeric",
                              sampler.result="data.frame",
                              marginal.prob="numeric",
                              sampled.res="list"))


bard.sampler.class<-function(bard.result,gamma,num_draws,sampler.result,marginal.prob,sampled.res,...)
{
    .bard.sampler.class(bard.result=bard.result,
                        gamma=gamma,
                        num_draws=num_draws,
                        sampler.result=sampler.result,
                        marginal.prob=marginal.prob,
                        sampled.res=sampled.res)
}



#' Detection of multivariate anomalous segments using BARD.
#'
#' Implements the BARD (Bayesian Abnormal Region Detector) procedure of Bardwell and Fearnhead (2017). 
#' BARD is a fully Bayesian inference procedure which is able to give measures of uncertainty about the 
#' number and location of anomalous regions. It uses negative binomial prior distributions on the lengths 
#' of anomalous and non-anomalous regions as well as a uniform prior for the means of anomalous regions. 
#' Inference is conducted by solving a set of recursions. To reduce computational and storage costs a resampling 
#' step is included.
#' 
#' @param x An n x p real matrix representing n observations of p variates.
#' @param p_N Hyper-parameter of the negative binomial distribution for the length of non-anomalous segments (probability of success). Defaults to \deqn{\frac{1}{n+1}.}
#' @param p_A Hyper-parameter of the negative binomial distribution for the length of anomalous segments (probability of success). Defaults to \deqn{\frac{5}{n}.}
#' @param k_N Hyper-parameter of the negative binomial distribution for the length of non-anomalous segments (size). Defaults to 1.
#' @param k_A Hyper-parameter of the negative binomial distribution for the length of anomalous segments (size). Defaults to \deqn{\frac{5p_A}{1- p_A}.}
#' @param pi_N Probability that an anomalous segment is followed by a non-anomalous segment. Defaults to 0.9.
#' @param paffected Proportion of the variates believed to be affected by any given anomalous segment. Defaults to 5\%. 
#' This parameter is relatively robust to being mis-specified and is studied empirically in Section 5.1 of \insertCite{bardwell2017;textual}{anomaly}.
#' @param lower The lower limit of the the prior uniform distribution for the mean of an anomalous segment \eqn{\mu}. Defaults to \deqn{2\sqrt{\frac{\log(n)}{n}}.}
#' @param upper The upper limit of the prior uniform distribution for the mean of an anomalous segment \eqn{\mu}. 
#' Defaults to the largest standardised value of x, i.e. \code{max(transform(x))}.
#' @param alpha Threshold used to control the resampling in the approximation of the posterior distribution at each time step. A sensible default is 1e-4.
#' Decreasing alpha increases the accuracy of the posterior distribution but also increases the computational complexity of the algorithm. 
#' @param h The step size in the numerical integration used to find the marginal likelihood. 
#' The quadrature points are located from \code{lower} to \code{upper} in steps of \code{h}. Defaults to 0.25. 
#' Decreasing this parameter increases the accuracy of the calculation for the marginal likelihood but increases computational complexity.   
#' @param transform A function used to transform the data prior to analysis. The default value is to scale the data using the median and the median absolute deviation.
#' 
#' @section Notes on default hyper-parameters:
#' This function gives certain default hyper-parameters for the two segment length distributions.
#' We chose these to be quite flexible for a range of problems. For non-anomalous segments a geometric distribution
#' was selected having an average segment length of \eqn{n} with the standard deviation being of the same order. 
#' For anomalous segments we chose parameters that gave an average length of 5 and a variance of \eqn{n}. 
#' These may not be suitable for all problems and the user is encouraged to tune these parameters. 
#' 
#' @return An instance of the S4 object of type \code{.bard.class} containing the data \code{x}, procedure parameter values, and the results.
#'
#' @references  \insertRef{bardwell2017}{anomaly}
#'
#' @seealso \code{\link{sampler}}
#'
#' @examples
#' 
#' library(anomaly)
#' set.seed(0)
#' sim.data<-simulate(n=500,p=50,mu=2,locations=c(100,200,300),
#'                    duration=6,proportions=c(0.04,0.06,0.08))
#' # run bard
#' bard.res<-bard(sim.data, alpha = 1e-3, h = 0.5)
#' sampler.res<-sampler(bard.res)
#' collective_anomalies(sampler.res)
#' \donttest{
#' plot(sampler.res,marginals=TRUE)
#' }
#' @export
bard<-function(x, p_N = 1/(nrow(x)+1), p_A = 5/nrow(x), k_N = 1, k_A = (5*p_A)/(1-p_A), pi_N = 0.9, paffected = 0.05, lower = 2*sqrt(log(nrow(x))/nrow(x)), upper = max(transform(x)), alpha=1e-4, h=0.25, transform=robustscale)
{
    # check the data
    x <- as.array(as.matrix(x))
    if(!is_array(x))
    {
        stop("cannot convert data to an array")
    }
    if(!all(is_not_na(x)))
    {
        stop("x contains NA values")
    }
    if(!all(is_not_null(x)))
    {
        stop("x contains NULL values")
    }
    if(any(is.infinite(x)))
    {
      stop("x contains Inf values")
    }
    if(!is_numeric(x))
    {
        stop("x must be of type numeric")
    }
    if(!is_function(transform))
    {
      stop("transform must be a function")
    }
    # transform the data
    x <- transform(x)
    
    # now convert the data to a list of vectors for marshalling to Rcpp
    # data<-Map(function(i) unlist(data[i,]),1:nrow(data))

    # check p_N
    if(!check.p(p_N))
    {
        stop("p_N must be in the interval (0,1)")
    }
    
    # check p_A
    if(!check.p(p_A))
    {
        stop("p_A must be in the interval (0,1)")
    }
    
    # check pi_N
    if(!check.p(pi_N))
    {
        stop("pi_N must be in the interval (0,1)")
    }
    
    # check alpha
    if(!check.p(alpha))
    {
      stop("alpha must be in the interval (0,1)")
    }
    
    # check paffected
    if(!check.p(paffected))
    {
        stop("paffected must be in the interval (0,1)")
    }
    
    # check k_A
    if(!check.k(k_A))
    {
      stop("k_A must be a positive real number")
    }
    
    # check k_N
    if(!check.k(k_N))
    {
      stop("k_N must be a positive real number")
    }
    
    # check lower
    if(!check.lu(lower))
    {
        stop("lower must be a real number")
    }

    # check upper
    if(!check.lu(upper))
    {
      stop("upper must be a real number")
    }
    
    # check relationaship between upper an lower
    if(lower > upper)
    {
        stop("value of upper should be greater than lower")
    }

    # check h
    if(!check.lu(h))
    {
      stop("h must be a real number")
    }
    
    # check relationship between h and upper and lower
    if(h > (upper - lower))
    {
        stop("h value too large : (h <= upper -lower)")
    }



#  CALL LB's master code !!!!!
# set up parameters
bardparams<-c(k_N, p_N, k_A, p_A, pi_N, paffected*dim(x)[2])
mu_seq<-c(seq(-upper, -lower, by=h), seq(lower, upper, by=h))
res<-Rbard(x, bardparams, mu_seq, alpha)

return(bard.class(data=x,
                  p_N=p_N,
                  p_A=p_A,
                  k_N=k_N,
                  k_A=k_A,
                  pi_N=pi_N,
                  alpha=alpha,
                  paffected=paffected,
                  lower=lower,
                  upper=upper,
                  h=h,
                  Rs=res)
           )




#    # set boost seed (used by c++) to sink with R random number generator
#    seed = as.integer(runif(1,0,.Machine$integer.max))

#    # dispatch    
#    Rs<-marshall_bard(data,p_N,p_A,k_N,k_A,pi_N,alpha,paffected,lower,upper,h,seed)

#    # put the data back into an array
#    data<-t(array(unlist(data),c(length(data[[1]]),length(data))))
    

#    return(bard.class(data=data,
#                      p_N=p_N,
#                      p_A=p_A,
#                      k_N=as.integer(k_N),
#                      k_A=as.integer(k_A),
#                      pi_N=pi_N,
#                      alpha=alpha,
#                      paffected=paffected,
#                      lower=lower,
#                      upper=upper,
#                      h=h,
#                      Rs=Rs)
#           )

    
}



#### helper functions for post processing

get.probvec <- function( filtering , params , num_draws=1000 )
{
  logprobs <- filtering$weights
  cpt.locs <- filtering$locations
  types <- filtering$type
  n <- length(filtering$weights)
  vec<-numeric(n)
  prob.states <- matrix(nrow=num_draws,ncol=n)
  
  for (i in 1:num_draws){
  
    f <- draw.from.post(logprobs, cpt.locs, types, params)
    segs <- f$draws
    states <- f$states
  
    # fill vec (hist of chpts)
    vec[f$draws] <- vec[f$draws] + 1
  
    # probs of being N or A
    prob.states[i,] <- get.states(segs,states,n)
    
  }

  return( apply(prob.states,2,sum)/num_draws )

}


# function to get a vector of states for a drawn value frin filtering posterior # 
# segs is vector of locations of segments
# states vector 
get.states <- function(segs, states, n){
  
  state.vec <- numeric(n)
  state.vec[ segs[1]:segs[2] ] <- states[1]
  if (length(segs) > 2){
    for (i in 2:( length(segs) - 1 ) ){
      state.vec[ ( segs[i]-1 ):segs[i+1] ] <- states[i]
    }
  }
  return(state.vec)
  
}


format_output = function(R){
  
  n = length(R)
  weights = vector("list", n) 
  location = vector("list", n)
  type = vector("list", n)
  for (i in 1:n){
    
    m = length(R[[i]])
    tmpweights = numeric(2*m)
    tmplocs = numeric(2*m)
    tmptypes = numeric(2*m)
    for (j in 1:m){
      triple = R[[i]][[j]]
      tmpweights[(2*j-1):(2*j)] = triple[1:2] 
      tmplocs[(2*j-1):(2*j)] = as.integer(triple[3])
      tmptypes[(2*j-1):(2*j)] = c(0,1) 
    }
    weights[[i]] = tmpweights[tmpweights > -Inf]
    location[[i]] = tmplocs[tmpweights > -Inf]
    type[[i]] = tmptypes[tmpweights > -Inf]
    
  }
  
  filtering = list("weights"=weights, "locations"=location, "type"=type)
  return(filtering)
  
}






#' Post processing of BARD results.
#'
#' Draw samples from the posterior distribution to give the locations of anomalous segments.
#'
#' @param bard_result An instance of the S4 class \code{.bard.class} containing a result returned by the \code{bard} function. 
#' @param gamma Parameter of loss function giving the cost of a false negative i.e. incorrectly allocating an anomalous point as being non-anomalous. 
#' For more details see Section 3.5 of \insertCite{bardwell2017;textual}{anomaly}.
#' @param num_draws Number of samples to draw from the posterior distribution. 
#' 
#' @return Returns an S4 class of type \code{bard.sampler.class}.  
#'
#' @references  \insertRef{bardwell2017}{anomaly}
#' @seealso \code{\link{bard}}
#'
#'
#' @examples
#' library(anomaly)
#' set.seed(0)
#' sim.data<-simulate(n=500,p=50,mu=2,locations=c(100,200,300),
#' duration=6,proportions=c(0.04,0.06,0.08))
#' # run bard
#' res<-bard(sim.data, alpha = 1e-3, h = 0.5)
#' # sample 
#' sampler(res)
#' 
#' @export
sampler<-function(bard_result, gamma = 1/3, num_draws = 1000)
{
   res<-Rsampler(res=bard_result@Rs,gamma=gamma,no.draws=num_draws)
   return(bard.sampler.class(bard_result,gamma,num_draws,res[[1]],res[[2]],res[[3]]))
}


#' @name plot-bard.sampler.class
#'
#' @docType methods
#'
#' @param marginals Logical value. If \code{marginals=TRUE} the plot will include visualisations of the marginal probablities of each time point being anomalous.
#' The defualt is \code{marginals=FALSE}. 
#' 
#' @rdname plot-methods
#'
#' @aliases plot,bard.sampler.class-method
#'
#' @export
setMethod("plot",signature=list("bard.sampler.class"),function(x,subset,variate_names,tile_plot,marginals=FALSE)
{   # from a plotting perspective the sampled BARD results are identical to those of PASS - cast type to pass.class
    # NULL out ggplot variables so as to pass CRAN checks
    variable<-prob<-value<-NULL
    tile.plot<-plot(pass.class(x@bard.result@data,x@sampler.result,0,0,0,0))
    tile.plot<-tile.plot+theme_bw()
    tile.plot<-tile.plot+theme(axis.ticks.y=element_blank())
    tile.plot<-tile.plot+theme(axis.text.y=element_blank(),axis.title=element_blank())
    tile.plot<-tile.plot+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    tile.plot<-tile.plot+xlab(label="t") # does not seem to shoe - need to find a fix

    if(marginals == FALSE)
    {
        return(tile.plot)
    }
    tile.plot<-tile.plot+theme(legend.position="none")
    df<-data.frame("t"=1:length(x@marginal.prob),"prob"=x@marginal.prob)
    marginal.prob.plot<-ggplot(df,aes(x=t,y=prob))
    marginal.prob.plot<-marginal.prob.plot+geom_line()
    marginal.prob.plot<-marginal.prob.plot+geom_hline(yintercept=1.0/(1.0+x@gamma),linetype="dashed",color="red")
    marginal.prob.plot<-marginal.prob.plot+theme_bw()
    marginal.prob.plot<-marginal.prob.plot+theme(axis.text.y=element_blank())
    marginal.prob.plot<-marginal.prob.plot+labs(x="t")
    marginal.prob.plot<-marginal.prob.plot+theme(strip.text.y=element_blank())


   gen.sample.plot<-function(object)
   {
       n<-nrow(object@bard.result@data) 
       df<-Reduce(cbind,
                  Map(function(locations)
                  {
                    x<-rep(0,n)
                    if (ncol(locations) > 0){
                      for(j in 1:ncol(locations))
                      {
                        s<-locations[1,j]
                        e<-locations[2,j]	
                        x[s:e]<-1
                      } 
                    }
                    return(x)
                  },
                  object@sampled.res)
       )
       n.df<-data.frame("n"=seq(1,nrow(df)))
       molten.data<-melt(cbind(n.df,df),id="n")
       out<-ggplot(molten.data, aes(n,variable))
       out<-out+geom_tile(aes(fill=value))
       #out<-out + theme(legend.position="none")
       #out<-out + theme(axis.text.y=element_blank())
       #out<-out + theme(axis.title.y=element_blank())
       #out<-out + theme(axis.line.y=element_blank())
       #out<-out + theme(axis.ticks.y=element_blank())
       #out<-out + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

       out<-out+theme_bw()
       out<-out+theme(axis.ticks.y=element_blank())
       out<-out+theme(axis.text.y=element_blank(),axis.title=element_blank())
       out<-out+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
       out<-out + theme(legend.position="none")
       
       return(out)
   }
   
    sample.plot<-gen.sample.plot(x)
    return( suppressWarnings(cowplot::plot_grid(tile.plot,sample.plot,marginal.prob.plot, align = "v", ncol = 1, rel_heights = c(1, 1,1))) )
})


#' @name show
#'
#' @docType methods
#'
#' @rdname show-methods
#'
#' @aliases show,bard.class-method
#'
#' @export
setMethod("show",signature=list("bard.class"),function(object)
{
    cat("BARD detecting changes in mean","\n",sep="")
    cat("observations = ",dim(object@data)[1],sep="")
    cat("\n",sep="")
    cat("variates = ",dim(object@data)[2],"\n",sep="")
    invisible()
})


#' @name summary
#'
#' @docType methods
#'
#' @rdname summary-methods
#'
#' @aliases summary,bard.sampler.class-method
#'
#' @export
setMethod("summary",signature=list("bard.sampler.class"),function(object,...)
{
    cat("BARD sampler detecting changes in mean","\n",sep="")
    cat("observations = ",dim(object@bard.result@data)[1],sep="")
    cat("\n",sep="")
    cat("variates = ",dim(object@bard.result@data)[2],"\n",sep="")
    cat("Collective anomalies detected : ",nrow(object@sampler.result),"\n",sep="")
    invisible()
})


#' @name show
#'
#' @docType methods
#'
#' @rdname show-methods
#'
#' @aliases show,bard.sampler.class-method
#'
#' @export
if(!isGeneric("collective_anomalies")) {setGeneric("collective_anomalies",function(object,...) {standardGeneric("collective_anomalies")})}
setMethod("show",signature=list("bard.sampler.class"),function(object)
{
    summary(object)
    ## LB changed - no anom doesnt show data frame with 0 columns and 0 rows
    if (nrow(object@sampler.result) > 0){
      print(object@sampler.result)
    }
    invisible()
})

#' @name collective_anomalies
#'
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,bard.sampler.class-method
#'
#' 
#' 
#' @export
setMethod("collective_anomalies",signature=list("bard.sampler.class"),function(object)
{
    return(object@sampler.result)
})


check.k = function(input){
  
  res = (length(input) == 1) && (is.numeric(input)) && (!is.nan(input)) && (!is.infinite(input)) && (input > 0)
  return(res)
  
}

check.p = function(input){
  
  res = (length(input) == 1) && (is.numeric(input)) && (!is.nan(input)) && (!is.infinite(input)) && (input > 0) && (input < 1)
  return(res)
  
}

check.lu = function(input){
  
  res = (length(input) == 1) && (is.numeric(input)) && (!is.nan(input))
  return(res)
  
}
