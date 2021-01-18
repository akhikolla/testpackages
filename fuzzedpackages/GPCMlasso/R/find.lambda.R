find.lambda <- function(cur.lambda, l.lambda, alpha.start, log_score_fun,
                        Q,q,I,n,response,design,designX, px, 
                        GHweights, GHnodes, acoefs, scale_fac, lambda2, cvalue, 
                        cores, weight, n_sigma, null_thresh, gradtol, iterlim, 
                        steptol){
  
  
  cur.low <- 0
  cur.up <- lambda.ratio <- Inf
  cur.lambda <- 0.5*cur.lambda

  while(lambda.ratio>0.15){
    suppressWarnings(
    m.opt <- try(nlm(log_score_fun, alpha.start,
                     Q = Q, q = q, I = I, n = n, Y = response, X = design, Z = designX, px = px,
                     GHweights = GHweights, GHnodes = GHnodes,
                     acoefs = acoefs, lambda = cur.lambda, scale_fac = scale_fac,
                     lambda2 = lambda2, cvalue = cvalue, cores = cores, weight = weight,  n_sigma = n_sigma,
                     gradtol = gradtol, iterlim = iterlim, check.analyticals = FALSE, 
                     steptol = steptol)))
  if(inherits(m.opt, "try-error")){
  stop("Determination of maximal tuning parameter failed. Please pre-specify
          the vector of tuning parameters manually via ctrl_GPCMlasso!")
    }

  coefs.cur <- m.opt$estimate%*%acoefs
  coefs.cur[abs(coefs.cur) < null_thresh] <-0

  if(sum(abs(coefs.cur))==0){
    cur.up <- cur.lambda
  }else{
    cur.low <- cur.lambda
  }
  if(is.finite(cur.up)){
    cur.lambda <- (cur.up-cur.low)*0.5+cur.low
    lambda.ratio <- (cur.up-cur.low)/cur.up
  }else{
    cur.lambda <- cur.lambda*2
  }

  ##end of while loop
  }
  return(cur.lambda)
}