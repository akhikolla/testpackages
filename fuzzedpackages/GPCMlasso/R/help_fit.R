loglikscoreDIFlasso2 <- function(alpha, Y, X, Z, Q, q, n, I,
                                 px, GHweights, GHnodes,
                                 acoefs, lambda, lambda2, cvalue, cores,
                                 weight, n_sigma, scale_fac){
  l <- loglikscoreDIFlasso(alpha, Y, X, Z, Q, q, n, I,
  px, GHweights, GHnodes,
  acoefs, lambda, lambda2, cvalue, cores,
  weight, n_sigma, scale_fac)
  
  ret <- l$objective
  attr(ret,"gradient") <- l$gradient
  ret
}

loglikscorePCMlasso2 <- function(alpha, Y, X, Z, Q, q, n, I,
                                 px, GHweights, GHnodes,
                                 acoefs, lambda, lambda2, cvalue, cores,
                                 weight, n_sigma, scale_fac){
  l <- loglikscorePCMlasso(alpha, Y, X, Z, Q, q, n, I,
                           px, GHweights, GHnodes,
                           acoefs, lambda, lambda2, cvalue, cores,
                           weight, n_sigma, scale_fac)
  
  ret <- l$objective
  attr(ret,"gradient") <- l$gradient
  ret
}

help_fit <- function(model, Y, l.lambda, start, loglik_fun, score_fun, log_score_fun, adaptive,
                        Q,q,I,n,m,response,design,designX, px,
                        GHweights, GHnodes, acoefs, lambda2, cvalue, n_sigma,
                        l.bound, trace, log.lambda, weight.penalties, scale_fac = scale_fac,
                        ada.lambda, lambda.min, ada.power, cores,
                        null_thresh, DSF, gradtol,iterlim, steptol, main.effects){


  ## get initial weight parameters   
  weight <- rep(1,ncol(acoefs))

  ## initialize starting values
  if(is.null(start)){
    if(trace){
      cat("Find start values ...\n")
    }
    alpha.start <- c(rep(0.1, px))

    if(!(model %in% c("RSM","GRSM"))){
    if(model=="GPCM"){
      m.ltm <- gpcm(Y, constraint = "gpcm")
      coef.ltm <- coef(m.ltm)
      if(is.matrix(coef.ltm)){
        sigma.start <- coef.ltm[,q[1]+1]
        delta.start <- c(t(coef.ltm[,-(q[1]+1)]))
      }else{
        sigma.start <- delta.start <- c()
        for(u in 1:I){
          sigma.start[u] <- coef.ltm[[u]][q[u]+1]
          delta.start <- c(delta.start, coef.ltm[[u]][-(q[u]+1)])
        }
      }
      
    }
    if(model=="PCM"){
      m.ltm <- gpcm(Y, constraint = "1PL")
      coef.ltm <- coef(m.ltm)
      if(is.matrix(coef.ltm)){
        sigma.start <- coef.ltm[1,q[1]+1]
        delta.start <- c(t(coef.ltm[,-(q[1]+1)]))
      }else{
        delta.start <- c()
        sigma.start <- coef.ltm[[1]][q[1]+1]
        for(u in 1:I){
          delta.start <- c(delta.start, coef.ltm[[u]][-(q[u]+1)])
        }
      }
    }
    if(model=="2PL"){
      m.ltm <- gpcm(Y, constraint = "gpcm")
      coef.ltm <- coef(m.ltm)
      sigma.start <- coef.ltm[,q[1]+1]
      delta.start <- c(t(coef.ltm[,-(q[1]+1)]))
    }
    if(model=="RM"){
      m.ltm <- gpcm(Y, constraint = "1PL")
      coef.ltm <- coef(m.ltm)
      sigma.start <- coef.ltm[1,q[1]+1]
      delta.start <- c(t(coef.ltm[,-(q[1]+1)]))
    }
    
    }else{

    if(model=="RSM"){
      m.mirt <- mirt(Y, 1, itemtype = 'rsm', verbose = FALSE)
      coefmethod <- selectMethod("coef", class(m.mirt))
      coef.mirt <- coefmethod(m.mirt, simplify = TRUE)
      sigma.start <- sqrt(coef.mirt$cov)
      alpha.mirt <- (coef.mirt$items)[1,2:(q[1]+1)]
      delta.mirt <- (coef.mirt$items)[,(q[1]+2)]
      delta.mirt <- -delta.mirt + alpha.mirt[1]
      alpha.mirt <- alpha.mirt - alpha.mirt[1]
      delta.start <- c(delta.mirt, alpha.mirt[-1])
    }
      if(model=="GRSM"){
        m.mirt <- mirt(Y, 1, itemtype = 'rsm', verbose = FALSE)
        coefmethod <- selectMethod("coef", class(m.mirt))
        coef.mirt <- coefmethod(m.mirt, simplify = TRUE)
        sigma.start <- rep(sqrt(coef.mirt$cov),n_sigma)
        alpha.mirt <- (coef.mirt$items)[1,2:(q[1]+1)]
        delta.mirt <- (coef.mirt$items)[,(q[1]+2)]
        delta.mirt <- -delta.mirt + alpha.mirt[1]
        alpha.mirt <- alpha.mirt - alpha.mirt[1]
        delta.start <- c(delta.mirt, alpha.mirt[-1])
      }
    }
    
      alpha.start <- c(delta.start,rep(0,ncol(designX)),abs(sigma.start))

      alpha.null <- alpha.start[rowSums(abs(acoefs))==0]
      p_null <- length(alpha.null)
      design_null <- matrix(0,0,0)
      if(main.effects & ncol(designX)>0){
        design_null <- designX[,1:m, drop = FALSE]
      }
      acoefs_null <- matrix(0,nrow=p_null,ncol=1)
      bound_null <- l.bound[rowSums(abs(acoefs))==0]

      loglik_NA <- TRUE

      while(loglik_NA){
        loglik_NA <- is.nan(log_score_fun(alpha.null,
                                Q = Q, q = q, I = I, n = n, Y = response, X = design,
                                Z = design_null, px = p_null,
                                GHweights = GHweights, GHnodes = GHnodes,
                                acoefs = acoefs_null, lambda = 0, scale_fac = scale_fac,
                                lambda2 = lambda2, cvalue = cvalue, cores = cores, weight = 1,  n_sigma = n_sigma))
        if(loglik_NA){
          alpha.null <- alpha.null*0.9
        }
              
      }

      m.opt <- try(nlminb(start = alpha.null, objective = loglik_fun, gradient = score_fun,
                    Q = Q, q = q, I = I, n = n, Y = response, X = design,
                    Z = design_null, px = p_null,
                    GHweights = GHweights, GHnodes = GHnodes,
                    acoefs = acoefs_null, lambda = 0, scale_fac = scale_fac,
                    lambda2 = lambda2, cvalue = cvalue, cores = cores, weight = 1,  n_sigma = n_sigma,
                    control=list(eval.max=500, iter.max=500, step.min=0.01),
                    lower=bound_null))

      alpha.start[rowSums(abs(acoefs))==0] <- m.opt$par
      alpha.start[rowSums(abs(acoefs))!=0] <- 1e-8

  }else{
    alpha.start <- start
  }

  
  ##  get new weights if necessary
  if(adaptive){
    if(trace){
      cat("Get adaptive weights ...", "\n")
    }
 

    m.opt <- try(nlminb(start = alpha.start, objective = loglik_fun, gradient = score_fun,
                        Q = Q, q = q, I = I, n = n, Y = response, X = design,
                        Z = designX, px = px,
                        GHweights = GHweights, GHnodes = GHnodes,
                        acoefs = acoefs, lambda = 0, scale_fac = scale_fac,
                        lambda2 = ada.lambda, cvalue = cvalue, cores = cores, weight = weight,  n_sigma = n_sigma,
                        control=list(eval.max=500, iter.max=500,step.min=0.01),
                        lower=l.bound))
    
    weight <- try(m.opt$par)

    if(inherits(m.opt,"try-error")){
      stop("Adaptive weights can not be calculated! Increase ada.lambda or set adaptive = FALSE!")
    }
    weight <- abs(t(acoefs) %*% weight)^ada.power
    if (any(weight == 0)) 
      weight[which(weight == 0)] <- 1e-4
    weight <- as.vector(1/weight)
  }

  if(weight.penalties & length(weight)!=(I*m) & m>0){
    weight.i <- (q+choose(q,2))/q
    weight.help <- c(rep(weight.i,q*m), rep(weight.i,choose(q,2)*m))
    weight.help <- weight.help/(prod(weight.help)^(1/length(weight.help)))
    weight <- weight/weight.help
  }
  
  #### find maximal lambda value and make grid
  if(!is.na(l.lambda)){
    
    if(trace){
      cat("Find maximal tuning parameter ...", "\n")
    }
    
  score <- score_fun(alpha.start, response, design, designX, Q,
                         q, n, I, px, GHweights, GHnodes, acoefs,
                         0, lambda2, cvalue,cores,  weight,
                         n_sigma, scale_fac)

  a <- abs(score/acoefs%*%weight)
  a[a==Inf] <- 0
  lambda.max <- max(a[rowSums(abs(acoefs))!=0])*1.1

    if(DSF){
    lambda.max <- find.lambda(lambda.max, l.lambda, alpha.start, log_score_fun,
                              Q,q,I,n,response,design,designX, px, 
                              GHweights, GHnodes, acoefs, scale_fac, lambda2, cvalue, 
                              cores, weight, n_sigma, null_thresh, gradtol, iterlim, 
                              steptol)
  }

if(log.lambda){
    lambda <- exp(seq(log(lambda.max+0.05*lambda.max), 
                      log(lambda.min+0.05*lambda.max), length = l.lambda))-0.05*lambda.max
    lambda[l.lambda] <- lambda.min
  }else{
    lambda <- seq(lambda.max,lambda.min,length=l.lambda)
  }
  }else{
    lambda <- NA
  }

  return(list(lambda = lambda, weight = weight, alpha.start = alpha.start))
}