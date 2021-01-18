find.lambda <- function(response, design, penalties, 
                        k, m, control, trace){
  
  if(trace){
    cat("Find maximal lambda...\n")
  }

  cur.low <- 0
  cur.up <- lambda.ratio <- Inf

  cur.lambda <- 100
  
  while(lambda.ratio>0.15){

    m.cur <- try(fit.BTLLasso(response, design, penalties, 
                    cur.lambda, k, m, control, trace = FALSE))


  # coefs.cur <- m.cur$coefs[rowSums(abs(penalties$acoefs))!=0]
  coefs.cur <-m.cur$coefs%*%penalties$acoefs
  coefs.cur[abs(coefs.cur) < 1/(10^control$precision)] <-0

  if(sum(abs(coefs.cur))==0){
    cur.up <- cur.lambda
  }else{
    cur.low <- cur.lambda
  }
  if(is.finite(cur.up)){
    cur.lambda <- (cur.up-cur.low)*0.7+cur.low
    lambda.ratio <- (cur.up-cur.low)/cur.up
  }else{
    cur.lambda <- cur.lambda*2
  }

  ##end of while loop
  }

  if(control$log.lambda){
    lambda <- exp(seq(log(cur.up+0.01*cur.up), 
                      log(control$lambda.min+0.01*cur.up), length = control$l.lambda))-0.01*cur.up
    lambda[control$l.lambda] <- control$lambda.min
  }else{
    lambda <- seq(cur.up,control$lambda.min,length=control$l.lambda)
  }
  
  
  return(lambda)
}