
###############################################
#main function
tfCox = function(dat, ord=0, alpha=1, lambda.seq=NULL, discrete=NULL,
                 n.lambda=30, lambda.min.ratio = 0.01, tol=1e-6,
                 niter=1000, stepSize=25, backtracking=0){
  call = match.call()
  fit = list()
  X = dat$X
  p = ncol(X)
  status = dat$status
  time = dat$time
  #check the inputs
  if (is.null(ncol(X))) X = matrix(X, ncol=1)
  if (any(time<=0)) stop("'time' should be positive")
  if (any(!(status %in% c(0,1)))) stop("'status' should be either 0 or 1")
  if (length(status) != nrow(X)) stop("'status' and 'X' do not have the same number of subjects")
  if (length(time) != nrow(X)) stop("'time' and 'X' do not have the same number of subjects")
  if (!is.null(lambda.seq) & any(lambda.seq<=0)) stop("Value for 'lambda' must be positive")
  if (ord<0 | round(ord)!=ord) stop("Value for 'ord' must be a positive integer")
  if (alpha<0 | alpha>1) stop("alpha must be > 0 and < 1 ")
  if (tol<=0) stop("Value for 'tol' must be positive")
  if (!is.null(discrete) & any(!(discrete %in% 1:p))) stop("'discrete' should be vectors of integer between 1 and p")
  
  if (is.null(lambda.seq)) {
    lam_max = max_lam(dat, ord, alpha, discrete)
    lambda.seq = exp(seq(log(lam_max), log(lam_max*lambda.min.ratio),length.out=n.lambda))
  }
  #make sure lambda.seq is decreasing
  lambda.seq = sort(lambda.seq, decreasing=TRUE)
  n.lambda = length(lambda.seq)
  #matrices and lists to store results
  fit$ord = ord
  fit$alpha = alpha
  fit$lambda.seq = lambda.seq
  fit$theta.list = list()
  fit$num.knots = numeric(n.lambda)
  fit$num.nonsparse = numeric(n.lambda)
  
  for (k in 1:n.lambda) {
    if (k == 1) {
      one.fit = tfCox.helper(dat, ord=ord, alpha=alpha, lambda=lambda.seq[k],
                             discrete=discrete, tol=tol, theta0=NULL, niter=niter,
                             stepSize=stepSize,backtracking = backtracking)
    }
    else {
      one.fit = tfCox.helper(dat, ord=ord, alpha=alpha, lambda=lambda.seq[k],
                             discrete=discrete, theta0=as.matrix(one.fit$theta),
                             tol=tol, niter=niter, stepSize=stepSize,backtracking = backtracking)
    }
    fit$theta.list[[k]] = one.fit$theta
    fit$num.knots[k] = one.fit$knots
    fit$num.nonsparse[k] = one.fit$non_zero
  }
  fit$dat = dat
  fit$call = call
  class(fit) = "tfCox"
  return(fit)
}

tfCox.helper = function(dat, ord=0, alpha=1, lambda, discrete=NULL, theta0=NULL, 
                        tol=1e-6, niter=1000, stepSize=25, backtracking=0) {
  
  X = as.matrix(dat$X)
  n = nrow(X)
  p = ncol(X)
  time = dat$time
  status = dat$status
  Perm = apply(X, 2, order) #permutation matrix
  Rank = apply(X, 2, rank) #rank matrix
  indx_time = order(time, -status)
  timeOrd = time[indx_time]
  tie = count_tie(timeOrd)
  thin = numeric(p)
  vec_xtol = numeric(p)
  for (j in 1:p) {
    x = X[Perm[,j],j]
    mindx = min(diff(x))
    x_tol = 1e-6*max(IQR(x),diff(range(x))/2)
    vec_xtol[j] = x_tol
    if (mindx < x_tol) thin[j] = 1
  }
  if (is.null(theta0)) theta0 = matrix(0, nrow=n, ncol=p)
  ndis = length(discrete)
  if (is.null(discrete)) discrete = -1 else discrete = discrete - 1
  
  fit = .Call('_tfCox_tfCox_onelambda', PACKAGE='tfCox', ord, alpha, lambda, 
              discrete, ndis, stepSize, X, theta0, time, 
              status, indx_time, tie, n, p, Perm, 
              Rank, thin, vec_xtol, tol, niter, backtracking)
  
  theta = fit$theta
  intercept = fit$intercept
  iter = fit$iter  
  #number of knots and proportion of non-zero features
  non_zero = 0
  knots = numeric(p)
  for (j in 1:p) {
    if (!(j %in% discrete)) {
      unique.x = X[!duplicated(X[,j]), j]
      unique.theta = theta[!duplicated(X[,j]), j]
      ordj = order(unique.x)
      order.x = unique.x[ordj]
      order.theta = unique.theta[ordj]
      if (length(order.x)>=ord+2) {
        indx = c(1,which(diff(order.x) >= vec_xtol[j]) + 1)
        Dj = get.Dk(x=order.x[indx], k=ord)
        knots[j] = sum( abs(Dj %*% order.theta[indx]) > 1e-3 )
      }
    }
    if (sum(abs(theta[,j]))>1e-2) non_zero = non_zero + 1/p
  }
  
  knots = sum(knots)
  return(list(ord=ord, theta=theta, knots=knots, non_zero=non_zero,
              lambda=lambda, alpha=alpha, intercept=intercept, iter=iter))
}


###############################################
#k-fold cross-validation
cv_tfCox = function(dat, ord=0, alpha=1, discrete=NULL, lambda.seq=NULL,
                    lambda.min.ratio=0.01, n.lambda=30, n.fold=5, seed=NULL,
                    tol=1e-6, niter=1000, stepSize=25, backtracking=0) {
  call = match.call()
  X = dat$X
  time = dat$time
  status = dat$status
  n = nrow(X)
  p = ncol(dat$X)
  #check the inputs
  if (any(time<=0)) stop("'time' should be positive")
  if (any(!(status %in% c(0,1)))) stop("'status' should be either 0 or 1")
  if (length(status) != nrow(X)) stop("'status' and 'X' do not have the same number of subjects")
  if (length(time) != nrow(X)) stop("'time' and 'X' do not have the same number of subjects")
  if (!is.null(lambda.seq) & any(lambda.seq<=0)) stop("Value for 'lambda' must be positive")
  if (ord<0 | round(ord)!=ord) stop("Value for 'ord' must be a positive integer")
  if (alpha<0 | alpha>1) stop("alpha must be > 0 and < 1 ")
  if (tol<=0) stop("Value for 'tol' must be positive")
  if (!is.null(discrete) & any(!(discrete %in% 1:p))) stop("'discrete' should be between 1:p")
  
  if (is.null(lambda.seq)) {
    lam_max = max_lam(dat, ord, alpha, discrete)
    lambda.seq = exp(seq(log(lam_max), log(lam_max*lambda.min.ratio),length.out=n.lambda))
  }
  
  n.lambda = length(lambda.seq)
  loss = matrix(NA, n.fold, n.lambda)
  if (!is.null(seed)) set.seed(seed)
  folds <- cut(sample(n),breaks=n.fold,labels=FALSE)
  
  testData <- list()
  trainData <- list()
  
  for(i in 1:n.fold){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData$time <- dat$time[testIndexes]
    testData$status <- dat$status[testIndexes]
    testData$X <- as.matrix(dat$X[testIndexes,])
    
    trainData$time <- dat$time[-testIndexes]
    trainData$status <- dat$status[-testIndexes]
    trainData$X <- as.matrix(dat$X[-testIndexes,])
    
    for (k in 1:n.lambda) {
      if (k == 1)
        one.fit = tfCox.helper(dat=trainData, ord=ord, alpha=alpha, lambda=lambda.seq[k],
                               discrete=discrete, theta0=NULL, tol=tol, niter=niter, 
                               stepSize=stepSize,backtracking=backtracking)
      else
        one.fit = tfCox.helper(dat=trainData, ord=ord, alpha=alpha, lambda=lambda.seq[k], 
                               discrete=discrete, theta0=as.matrix(one.fit$theta), 
                               tol=tol, niter=niter, stepSize=stepSize,
                               backtracking=backtracking)
      
      theta_hat = predict.tfCox.helper(ord, trainData$X, testData$X, one.fit$theta)
      loss[i,k] = negloglik(testData, theta_hat)
    }
  }
  
  mean.cv.error <- apply(loss,2,mean)
  se.cv.error <- apply(loss,2,sd)/sqrt(n.fold)
  indx = which.min(mean.cv.error)
  result = list(best.lambda = lambda.seq[indx] ,lambda.seq=lambda.seq,
                mean.cv.error=mean.cv.error, se.cv.error=se.cv.error,
                ord=ord, alpha=alpha, n.fold=n.fold, call=call)
  class(result) = "cv_tfCox"
  return(result)
}

###############################################
# count the tie events for ordered time
count_tie = function(t) {
  n=length(t)
  tie = 0
  count = 0
  cur = 1
  i = 2
  while (i <= n ) {
    if (t[i] == t[cur]) {
      tie[i] = 0
      count = count + 1
    }
    else {
      tie[cur] = count
      cur = i
      count = 0
    }
    if (i == n) tie[cur] = count
    i=i+1
  }
  # no ties: 0
  # ties: 1,2,...
  return(tie)
}


#calculate gradient
gradient = function(dat, theta) {
  X = dat$X
  n = nrow(X)
  p = ncol(X)
  time = dat$time
  status = dat$status
  indx = order(time, -status)
  timeOrd = time[indx]
  statusOrd = status[indx]
  tie = count_tie(timeOrd)
  firstObs = n
  for (u in 1:n) {
    if (statusOrd[u] == 1){
      firstObs = u
      break
    }
  }
  theta = as.matrix(theta[indx,])
  grad = rep(0,n)
  exp.linear = sapply(1:n, function(k) exp(sum(theta[k,])))
  sum.riskset = numeric(n)
  for (u in 1:n) sum.riskset[u] = sum(exp.linear[u:n])
  
  for (u in 1:n) {
    i = firstObs
    while (i <= u) {
      if (statusOrd[i]) {
        prob_ui = exp.linear[u]/sum.riskset[i]
        grad[u] = grad[u] + sum(statusOrd[i:(i+tie[i])]) * prob_ui
      }
      if (tie[i]) i = i + tie[i]
      i=i+1
    }
    grad[u] = grad[u]/n
    if (statusOrd[u]) grad[u] = grad[u] - 1/n
  }
  grad[indx] = grad
  return(grad=grad)
}


#calculate the loss
negloglik = function(dat, theta) {
  X = dat$X
  n = nrow(X)
  p = ncol(X)
  time = dat$time
  status = dat$status
  indx = order(time, -status)
  time = time[indx]
  status = status[indx]
  theta = as.matrix(theta[indx,])
  tie = count_tie(time)
  s = 0
  linear_term = numeric(n)
  for (i in 1:n) {
    linear_term[i] = sum(theta[i,])
  }
  
  i = 1
  while (i <= n) {
    if (status[i]) {
      s = s + sum(linear_term[i:(i+tie[i])]) - sum(status[i:(i+tie[i])]) *log(sum(exp(linear_term[i:n])))
    }
    if (tie[i]) i = i + tie[i]
    i = i+1
  }
  s = s*(-1/n)
  return(s)
}



###############################################
#predict theta from "tfCox" object
predict.tfCox = function(object, newX, which.lambda=1, ...) {
  fit = object
  if (class(fit)!="tfCox") stop("Provide 'fit' of class 'tfCox'")
  dat = fit$dat
  X = dat$X
  ord = fit$ord
  theta_fit = fit$theta.list[[which.lambda]]
  return(predict.tfCox.helper(ord, X, newX, theta_fit))
}

predict.tfCox.helper = function(ord, X, newX, theta_fit) {
  p = ncol(X)
  theta_pred = matrix(0, nrow(newX), p)
  for (i in 1:p) {
    ordi = order(as.vector(X[,i]))
    if (ord == 0)
      theta_pred[,i] = approx(X[ordi,i], theta_fit[ordi,i], xout=newX[,i],
                              method="constant",rule=2)$y
    else
      theta_pred[,i] = approx(X[ordi,i], theta_fit[ordi,i], xout=newX[,i],
                              method="linear",rule=2)$y
  }
  return(theta_pred)
}

###############################################
#simulate data
gen_theta = function(x, scenario, feature) {
  option = scenario*10+feature
  
  #scenario 1: all piecewise constant f_j
  if (option==11) {theta_j = 3 * (x<(-1)) + 2 + -7 * (x>0.5); theta_j = (theta_j - 0.1)/4.3232149}
  if (option==12) {theta_j = 12 * (x<(-0.2)) + -5 + 7 * (x>1.1); theta_j = (theta_j - 2.48)/4.8999837}
  if (option==13) {theta_j = 3 + -6 * (x<(-1.7) | x>0.8); theta_j = theta_j/3.0000150}
  if (option==14) {theta_j = -5 + 6 * (x>(-.7)) + 3 * (x>1.6); theta_j = (theta_j +.62)/3.4577044}
  
  #scenario 2: smooth f_j - SPAM paper functions
  if (option==21) {theta_j = -sin(1.5*x); theta_j = theta_j/0.6614151}
  if (option==22) {theta_j = x^3 + 1.5*(x-0.5)^2; theta_j = (theta_j - 3.500063)/4.8930086}
  if (option==23) {theta_j = -pnorm(x, mean=0.5, sd=0.8); theta_j = (theta_j + .4003183)/0.3874514}
  if (option==24) {theta_j = sin(exp(-0.5*x)); theta_j = (theta_j - .6195458)/0.2985293}
  
  #scenario 3: piecewise linear f_j
  if (option==31) {theta_j = ifelse(x < -1.25, 3*x+3.75, ifelse(x > 1.25, 2.5*x-3.125, 0)) ; theta_j = (theta_j+0.078125)/1.406973}
  if (option==32) {theta_j = ifelse(x < -0.5, -2*x-1, ifelse(x > 1,-2*x+5.75,2.5*x+1.25)); theta_j = (theta_j-2.0375)/1.063921}
  if (option==33) {theta_j = ifelse(x < -1, 2.5*x+3.5, ifelse(x > 1, 2.5*x-5.5, -2*x-1)); theta_j = (theta_j+1)/1.116169}
  if (option==34) {theta_j = ifelse(x < 0, 2*x, -2*x); theta_j = (theta_j+2.5)/1.443376}
  
  #scenario 4: 'large n' with functions with large constant functions
  if (option==41) {theta_j = -5 * (x<0) + (10*(x)^2 - 5) * (x>=0 & x<1) + (sin((x-1)*20) + 5) * (x>=1); theta_j = (theta_j + 1.324793)/4.562398}
  if (option==42) {theta_j = 2 + 3*(cos((x+.75)*2*pi) - 1) * (x>=-.75 & x<(-.25)) + (-4 + (cos((x+.75)*2*pi) - 1)) * (x>=-.25 & x<=0.25) + -4*(x>.25); theta_j = (theta_j + 0.59994)/2.083371}
  if (option==43) {theta_j = (x<(-1)) * 3.125 + (x>1) * -3.125 + (x>=-1 & x<(-0.5))*(-50*(x+1)^5 + 3.125) + (x>=-0.5 & x<.5)*-50*(x)^5 + (x>=0.5 & x<=1)*(-50*(x-1)^5 + -3.125); theta_j = (theta_j)/2.752587}
  if (option==44) {theta_j = (x<(-1)) * ((cos((x+1+pi)*10-9*pi)+1)*5) + (x>1.5)* (cos((x-1.5)*20)-1); theta_j = (theta_j - 1.244387)/3.03563}
  
  return(theta_j)
}


sim_dat <- function(n, zerof=0, scenario=1, scale=1, shape=1, censoring.rate=0.01,
                    n.discrete=0) {
  if (!scenario %in% c(1,2,3,4)) stop("scenario must be 1, 2, 3, or 4")
  if (n.discrete<0) stop("n.discrete must be a non-negative integer")
  p = 4 + n.discrete + zerof
  X = matrix(runif(n = n * p, min = -2.5, max = 2.5), nrow = n, ncol = p)
  theta = matrix(0,n,p)
  for (i in 1:4) {
    theta[,i] = gen_theta(X[,i], scenario, i)
    theta[,i] = theta[,i] - mean(theta[,i])
  }
  if (n.discrete) {
    for (j in 5:(5+n.discrete-1)) {
      X[,j] = (X[,j]>0)*1
      theta[X[,j]==1,j] = 0.5
      theta[X[,j]==0,j] = -0.5
      theta[,j] = theta[,j] - mean(theta[,j])
    }
  }
  
  v = runif(n)
  Tlat = - (log(v) / (scale * exp(apply(theta, 1, sum) )))^(1 / shape)
  C = rexp(n, rate=censoring.rate)
  time = pmin(Tlat, C)
  time = as.vector(time)
  status = as.numeric(Tlat <= C)
  dat = list(time=time, status=status, X=X, true_theta=theta)
  return(dat)
}


plot.sim_dat = function(x, which.predictor = NULL, n.plot = 4, ...) {
  
  dat = x
  X = dat$X
  p = ncol(X)
  true_theta = dat$true_theta
  if (!is.null(which.predictor) & any(which.predictor > p) | any(which.predictor <=0)) stop("which.predictor should be integer vector from 1 to p")
  
  k = 0
  if (!is.null(which.predictor)) {
    for (i in which.predictor) {
      ordi = order(as.vector(X[,i]))
      plot(X[ordi,i], true_theta[ordi,i], xlab=paste("predictor",i),
           ylab="true theta" )
    }
  }
  else{
    for (i in 1:p) {
      ordi = order(as.vector(X[,i]))
      plot(X[ordi,i], true_theta[ordi,i], xlab=paste("predictor",i),
           ylab="true theta" )
      k = k+1
      if (k >= n.plot) break
    }
  }
}

##############################################
#calculate the maximxum lambda

#function to generate D^{(1)} with dim. (n-1) x n
#note that get.D1(n=length(x)) gives the same output as get.Dk(x=x, k=0)
get.D1 = function(n) {
  cbind(-diag(n-1), 0) + cbind(0, diag(n-1))
}

#function to generate D^{(x, k+1)} where k is order of trend filtering
get.Dk = function(x, k) {
  if (k==0) {
    get.D1(n = length(x))
  } else {
    get.D1(n = length(x)-k) %*% get.Diag(x = x, k = k) %*% get.Dk(x = x, k = k-1)
  }
}

#function to generate diagonal matrix that is part of D^{(x,k+1)}
get.Diag = function(x, k) {
  n = length(x)
  diag(x = k/(x[(k+1):n]-x[1:(n-k)]))
}

#calculate the maximum lambda
max_lam = function(dat, ord = 0, alpha = 1, discrete=NULL) {
  X = dat$X
  n = nrow(X)
  p = ncol(X)
  Perm = apply(X, 2, order)
  Rank = apply(X, 2, rank)
  
  if (is.null(discrete))
    penalized = 1:p
  else
    penalized = (1:p)[-discrete]
  
  if (alpha == 1) {
    #calculate the solution when tuning parameter is large
    theta = matrix(0, n, p)
    
    if (ord == 0) {
      for (j in discrete) {
        uni_fit = survival::coxph(survival::Surv(dat$time, dat$status) ~ as.factor(dat$X[,j]))
        theta[,j] = uni_fit$linear.predictors
      }
    }
    
    if (ord > 0) {
        for (j in discrete) {
          uni_fit = survival::coxph(survival::Surv(dat$time, dat$status) ~ as.factor(dat$X[,j]))
          theta[,j] = uni_fit$linear.predictors
        }
        for (j in penalized) {
          uni_fit = survival::coxph(survival::Surv(dat$time, dat$status) ~ I(dat$X[,j]^ord))
          theta[,j] = uni_fit$linear.predictors
        }
      
    }
    ############################
    
    grad1 = gradient(dat, theta)
    temp = numeric(p)
    for (j in 1:p) {
      if (!(j %in% discrete)) {
        order.x = X[Perm[,j],j]
        order.grad = grad1[Perm[,j]]
        indx = c(1,which(diff(order.x) >= 1e-6) + 1)
        Dj = get.Dk(x=order.x[indx], k=ord)
        coef = lm(order.grad[indx] ~ t(Dj))$coef
        temp[j] = max(abs(coef))
      }
    }
    max1 = max(temp)
    return(max1)
  }
  
  else if (alpha == 0) {
    grad2 = gradient(dat, matrix(0, n, p))
    max2 = sqrt(sum(grad2^2))
    return(max2)
  }
  
  else {
    grad3 = gradient(dat, matrix(0, n, p))
    max3 = sqrt(sum(grad3^2))/(1-alpha)
    return(max3)
  }
}

###############################################
#plot the fitted functions
plot.tfCox = function(x, which.lambda=1, which.predictor = NULL, n.plot = 4, ...){
  
  fit = x 
  dat = fit$dat
  X = dat$X
  p = ncol(X)
  ord = fit$ord
  if (!(which.lambda %in% 1:length(fit$lambda.seq))) stop("Provide a valid 'index'")
  if (!is.null(which.predictor) & any(which.predictor > p) | any(which.predictor <=0)) stop("which.predictor should be integer vector from 1 to p")
  
  theta = fit$theta.list[[which.lambda]]
  k = 0
  if (!is.null(which.predictor)) {
    for (i in which.predictor) {
      ordi = order(as.vector(X[,i]))
      if (ord == 0)
        plot(X[ordi,i], theta[ordi,i], type="s", xlab=paste("predictor",i),
             ylab="Estimated log HR" )
      else
        plot(X[ordi,i], theta[ordi,i], type="l", xlab=paste("predictor",i),
             ylab="Estimated log HR" )
    }
  }
  else{
    for (i in 1:p) {
      if (sum(theta[,i]^2) > 0) {
        ordi = order(as.vector(X[,i]))
        if (ord == 0)
          plot(X[ordi,i], theta[ordi,i], type="s", xlab=paste("predictor",i),
               ylab="Estimated log HR" )
        else
          plot(X[ordi,i], theta[ordi,i], type="l", xlab=paste("predictor",i),
               ylab="Estimated log HR" )
        k = k+1
      }
      if (k >= n.plot) break
    }
    if (k == 0) cat("no plot is shown because all thetas are zero")
  }
}

plot.cv_tfCox = function(x, showSE=F, ...) {
  cv = x
  min.y = min(cv$mean.cv.error - cv$se.cv.error)
  max.y = max(cv$mean.cv.error + cv$se.cv.error)
  plot(x=1,type="n",xlim=c(min(log(cv$lambda.seq)),max(log(cv$lambda.seq))),
       ylim=c(min.y,max.y),ylab="Cross-validation Error",xlab="log(Lambda)")
  points(log(cv$lambda.seq), cv$mean.cv.error, cex=1.5, pch=16)
  if (showSE==T) arrows(x0=log(cv$lambda.seq),x1=log(cv$lambda.seq),
                        y0=cv$mean.cv.error-cv$se.cv.error,
                        y1=cv$mean.cv.error+cv$se.cv.error,
                        angle=90,length=0.05,code=3,lwd=1.6)
  abline(v=log(cv$best.lambda),lty=2)
}

################################################
#summarize the fit of tfCox
summary.tfCox = function(object, ...) {
  fit = object
  if (class(fit)!="tfCox") stop("Provide 'fit' of class 'tfCox'")
  cat("Call: \n")
  print(fit$call)
  args = match.call()
  
  cat("\nThe tfCox fit corresponds to:\n")
  cat(paste(fit$ord,"th order trend filtering", sep=""),"\n")
  cat(paste("alpha = ",fit$alpha,sep=""), "\n")
  cat("\nThe number of knots and the percent sparsity for tuning parameter lambda:\n")
  res = matrix(c(fit$num.knots, round(100*(1-fit$num.nonsparse),2)),
               nrow=length(fit$lambda.seq), ncol=2)
  rownames(res) = round(fit$lambda.seq,3)
  colnames(res) = c("knots", "sparsity")
  print(res)
}

summary.cv_tfCox = function(object, ...) {
  cv = object
  if (class(cv) != "cv_tfCox") stop("Provide 'cv' of class 'cv_tfCox'")
  cat("Call: \n")
  print(cv$call)
  args = match.call()
  
  cat("\nThe cv_tfCox fit correspond to:\n\n")
  cat(paste(cv$ord,"th order trend filtering", sep=""),"\n")
  cat("lambda: "); cat(round(cv$lambda.seq,3))
  cat("\n\nalpha: "); cat(cv$alpha)
  
  cat(paste("\n\nCross-validation with K=",cv$n.fold," folds was used to choose lambda.",sep=""))
  cat("\nLambda was chosen to be the value with the minimum CV error.\n")
  cat(paste("\nThe chosen lambda was ",round(cv$best.lambda,3),".\n",sep=""))
  cat(paste("\nUse 'plot(",args$cv,")' to plot CV error curve.\n",sep=""))
}

###############################################
#choose the tuning parameter lambda using training and testing dataset
tfCox_choose_lambda = function(dat, test_dat, ord = 0, alpha=1, discrete=NULL,
          lam_seq = NULL, nlambda=30, c=NULL, tol=1e-6, niter=1000,
          stepSize=25, backtracking=0) {
  
  time = dat$time
  status = dat$status
  X = dat$X
  n = nrow(X)
  p = ncol(X)
  if (is.null(c)) {
    c = 0.01
  }
  if (is.null(lam_seq)) {
    lam_max = max_lam(dat, ord = ord, alpha, discrete)
    lam_seq = exp(seq(log(lam_max), log(lam_max*c),length.out=nlambda))
  }
  #make sure lambda.seq is decreasing
  lam_seq = sort(lam_seq, decreasing=TRUE)
  nlam = length(lam_seq)
  loss = numeric(nlam)
  knots = numeric(nlam)
  paramfit = numeric(nlam)
  best_theta = matrix(0, n, p)
  best_loss = Inf
  
  for (k in 1:nlam) {
    if (k == 1) {
      res = tfCox.helper(dat=dat, ord = ord, alpha=alpha, lambda=lam_seq[k],
            discrete=discrete, tol=tol, theta0=NULL, niter=niter, stepSize=stepSize,
            backtracking = backtracking )
    }
    else {
      res = tfCox.helper(dat=dat, ord = ord, alpha=alpha, lambda=lam_seq[k],
            discrete=discrete, tol=tol, theta0=res$theta, niter=niter, stepSize=stepSize,
            backtracking = backtracking)
    }
    knots[k] = res$knots
    theta_hat = predict.tfCox.helper(ord, X, test_dat$X, res$theta)
    loss[k] = negloglik(test_dat, theta_hat)
    paramfit[k] = sum((test_dat$true_theta - theta_hat)^2)

    if (best_loss > loss[k]) {
      best_lambda = lam_seq[k]
      best_knots = res$knots
      best_loss = loss[k]
      best_paramfit = paramfit[k]
      best_theta = res$theta
      best_nonzero = res$non_zero
    }
  }
  
  result = list(dat=dat, ord=ord, lam_seq=lam_seq, loss=loss, knots=knots,
                paramfit=paramfit, best_lambda=best_lambda, best_knots=best_knots,
                best_paramfit=best_paramfit, best_theta=best_theta,
                best_nonzero = best_nonzero)
  return(result)
}

#predict from the best lambda fit
predict_best_lambda = function(cv, newX) {
  dat = cv$dat
  X = dat$X
  p = ncol(X)
  ord = cv$ord
  theta_fit = cv$best_theta
  theta_pred = matrix(0, nrow(newX), p)
  
  for (i in 1:p) {
    ordi = order(as.vector(X[,i]))
    if (ord == 0)
      theta_pred[,i] = approx(X[ordi,i], theta_fit[ordi,i], xout=newX[,i], method="constant",rule=2)$y
    else if (ord == 1)
      theta_pred[,i] = approx(X[ordi,i], theta_fit[ordi,i], xout=newX[,i], method="linear",rule=2)$y
  }
  return(theta_hat = theta_pred)
}

