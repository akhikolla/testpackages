#' @title
#' robust cross-validation
#'
#' @description compute rocv
#'
#
#'
#'
#'
#'
#' @export
cv.gam <- function(X, Y, init.mode=c("sLTS","RLARS","RANSAC"),lambda.mode="lambda0",
                   lmax=1,lmin=0.05, nlambda=50, fold=10,
                   ncores=1, gam=0.1, gam0=0.5, intercept="TRUE", alpha=1,
                   ini.subsamp=0.2, ini.cand=1000,
                   alpha.LTS=0.75, nlambda.LTS=40){


  #########################################
  # detect the data type and tuning parameters value
  #########################################
  if(is.matrix(X)=="FALSE" || is.matrix(Y)=="FALSE" ){
    stop(message="The data type of X and Y is only Matrix.")
  }
  if(gam <0 || gam0 <0   ){
    stop(message="Robust tuning parameters gam and gam0 are only positive.")
  }
  fold <- as.integer(fold)
  if( fold<1 || fold > nrow(X)  ){
    stop(message="The number of folds for RoCV is more than 1 and less than sample size.")
  }
  ini.cand <- as.integer(ini.cand)
  if(ini.cand <2  ){
    stop(message="The number of candidates of initial points is more than 2.")
  }
  nlambda <- as.integer(nlambda)
  if(nlambda < 1 ){
    stop(message="The number of grids for sparse tuning parameter is more than 1.")
  }
  nlambda.LTS <- as.integer(nlambda.LTS)
  if(nlambda < 1 ){
    stop(message="The number of grids for sLTS is more than 1.")
  }
  if(lmin < 0 || lmax <0){
    stop(message="lmin and lmax are positive.")
  }
  ncores <- as.integer(ncores)
  if(ncores < 1){
    stop(message="The number of cpu cores for parallel computing is more than 1.")
  }
  if( ini.subsamp <0 || ini.subsamp >1){
    stop(message="The proportion of subsample is 0 < ini.subsamp < 1")
  }

  #########################################
  # intercept
  #########################################

  if(intercept=="TRUE"){
    inter <-1
  }else {
    inter <-0
  }

  #########################################
  # choose a initial point
  #########################################
  init <- match.arg(init.mode)

  if(init=="RANSAC"){
    ransac <- gam.ini(X, Y, ini.subsamp, ncores, ini.cand, alpha, gam0)

    beta0.init <- inter*(ransac$beta0)
    beta.init <- as.matrix(ransac$beta)
    sigma.init <- ransac$sigma
    sigma.cv <- sigma.init
  }else if(init=="sLTS"){
    frac<- seq(0,lambda0(X,Y) ,by=(1.0/nlambda.LTS)*lambda0(X,Y))
    sLTS <- sparseLTS(Y~X,lambda=frac[-1],mode="lambda",alpha=alpha.LTS,ncores=ncores)

    beta0.init <- inter*(coef(sLTS)[1])
    beta.init <- as.matrix(coef(sLTS)[-1])
    sigma.init <- getScale(sLTS)
    sigma.cv <- sigma.init
  }else if(init =="RLARS"){
    RLARS <- rlars(Y~X)

    beta0.init <- inter*(coef(RLARS)[1])
    beta.init <- as.matrix(coef(RLARS)[-1])
    sigma.init <- getScale(RLARS)
    sigma.cv <- sigma.init
  }

  ########################################
  # regularization grids
  ########################################
  lmode <- match.arg(lambda.mode)

  if(lmode=="lambda0"){
      tmp <- lambda.ini(X, Y, beta0.init, beta.init, sigma.init ,gam)
      ratio <- 0.05
      lambda <- exp(seq(  log(ratio*tmp)  ,log(tmp)  ,length=nlambda))
  }else{
    lambda <- exp(seq( log(lmin)  ,log(lmax)  ,length=nlambda))
  }

  #######################################
  # caluculate Rocv
  #######################################

  id.cv <- rep(1:fold,length=nrow(X))
  list.cv <- 1:fold

  res.cvtmp <- rep(0.0, nlambda)
  res.cv <- rep(0.0, nlambda)

  for (cv.fold in 1:fold){
    X.tra <-  subset(X, id.cv %in% list.cv[-cv.fold])
    Y.tra <-  subset(Y, id.cv %in% list.cv[-cv.fold])

    X.test <- subset(X, id.cv %in% c(cv.fold))
    Y.test <- subset(Y, id.cv %in% c(cv.fold))


    for(k in 1:nlambda){
      res <- gam_reg(X.tra, Y.tra, beta.init, beta0.init, sigma.init, lambda[k], gam, glmnet, inter, alpha)
      res.cv[k] <- sum(dnorm(Y.test, res$beta0*(numeric(length(Y.test))+1)+X.test%*%res$beta, sigma.cv)^gam0)
    }
    res.cvtmp <- res.cvtmp + res.cv
  }

  tmp2 <- (-1/gam0)*log(res.cvtmp)
  tmp2 <- sort.list(tmp2)[1]

  res.reggam <- list()
  for(k in 1:nlambda){
    res <- gam_reg(X, Y, beta.init, beta0.init,sigma.init,  lambda[k], gam, glmnet, inter, alpha)
    res.reggam[[k]] <- list(beta0=res$beta0, beta=res$beta, sigma=res$sigma)
  }

  cat("\n")
  cat("init.mode is" , init, "\n")
  cat("The number of lambda is" ,nlambda ,"\n")
  cat("gamma is", gam, "gamma0 is",gam0 ,"\n")

  return(list(lambda=lambda, fit=res.reggam, Rocv= res.reggam[[tmp2]] ) )
}


gam.ini <- function(X, Y, ini.subsamp, cl.num, ini.cand, reg.alpha, gam0){

  cl.ini <- makeCluster(cl.num)
  registerDoParallel(cl.ini)
  on.exit(stopCluster(cl.ini))
  ini.sam <- floor(nrow(X)*ini.subsamp)

  p <- ncol(X)
  n <- nrow(X)
  tmp <- foreach(k =1:ini.cand, .combine=rbind, .packages='glmnet')%dopar%{

      t <- sample(n,ini.sam)
      x <- X[t,]
      y <- Y[t,]

      x_t <- X[setdiff(1:n,t),]
      y_t <- Y[setdiff(1:n,t),]

      res <- glmnet(x,y,family="gaussian",alpha=reg.alpha)
      sigma <- apply( abs(matrix(y,nrow =nrow(x),ncol=length(res$lambda))- matrix(res$a0,nrow =nrow(x),ncol=length(res$lambda),byrow=TRUE)-x%*%res$beta ),2, median )
      sigma <- 1.4826*sigma
      sigma <- t(matrix(sigma,ncol =nrow(x_t),nrow=length(res$lambda)))

      a <- matrix(y_t,nrow =nrow(x_t),ncol=length(res$lambda))- matrix(res$a0,nrow =nrow(x_t),ncol=length(res$lambda),byrow=TRUE)-x_t%*%res$beta

      tmp1 <- (-1/gam0)*log ( (1/nrow(x_t))*colSums((exp(-a^2/(2*sigma^2))/sqrt((2*pi*sigma^2)) )^gam0))-(gam0/(2*(1+gam0)))*log(sigma[1,]^2)
      tmp2 <-  which.min(tmp1)

      return(  c(tmp1[tmp2],res$a0[tmp2],res$beta[,tmp2], sigma[1,tmp2]) )
    }

  ini.ind <-  which.min(tmp[,1])
  return( list(beta0=unname(tmp[ini.ind, 2]), beta=unname(tmp[ini.ind, 3:(p+2)]), sigma=unname(tmp[ini.ind, (p+3)])) )
}




lambda.ini <- function(X, Y, beta0, beta, sigma, gam){

  tmp4 <- (numeric(nrow(X))+1)
  tmp1 <- beta0*tmp4
  tmp2 <- X%*%beta

  alpha <- exp(-gam*(Y-tmp1-tmp2)^2/(2*sigma^2) )
  alpha <- alpha/sum(alpha)

  tmp1 <- sum(alpha*(Y - tmp2))*tmp4

  Y.lasso <- sqrt(alpha)*(Y-tmp1)
  X.lasso <- diag(c(sqrt(alpha)))%*%X

  return (lambda0(X.lasso/sigma, Y.lasso/sigma, intercept=FALSE))
}


