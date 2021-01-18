## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github(repo="sentian/BOSSreg", subdir="r-package")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages(repo="BOSSreg", repos = "http://cran.us.r-project.org")

## ------------------------------------------------------------------------
n = 200 # Number of observations
p = 14 # Number of predictors
p0 = 6 # Number of active predictors (beta_j=0)
rho = 0.9 # Correlation between predictors
nrep = 1000 # Number of replications of y to be generated
SNR = 7 # Signal-to-noise ratio
seed = 65 # The seed for reproducibility

## ------------------------------------------------------------------------
library(MASS)
# Function to generate the data
# Columns of X have mean 0 and norm 1, y has mean 0
simu.data <- function(n, p, p0, rho, nrep, SNR, seed){
  # True beta
  beta = rep(0,p)
  beta = c(rep(c(1,-1),p0/2), rep(0,p-p0))
  names(beta) = paste0('X', seq(1,p))
  # Covariance matrix
  covmatrix  = matrix(0,nrow=p,ncol=p)
  diag(covmatrix) = 1
  for(i in 1:(p0/2)){
    covmatrix[2*i-1,2*i] = covmatrix[2*i,2*i-1] = rho
  }
  # Generate the predictors given the correlation structure
  set.seed(seed)
  x = mvrnorm(n,mu=rep(0,p),Sigma=covmatrix)
  x = scale(x,center=TRUE,scale=FALSE)
  colnorm = apply(x,2,function(m){sqrt(sum(m^2))}) 
  x = scale(x,center=FALSE,scale=colnorm) # standardization
  # Sigma calculated based on SNR
  sd = sqrt(t(beta/colnorm)%*%covmatrix%*%(beta/colnorm) / SNR)
  mu = x%*%beta 
  # Generate replications of y by fixing X
  y = matrix(rep(mu,each=nrep),ncol=nrep,byrow=TRUE) +      
    scale(matrix(rnorm(n*nrep,mean=0,sd=sd),nrow=n,ncol=nrep),center=TRUE,scale=FALSE)
  return(list(x=x, y=y, beta=beta, sigma=sd))
}
dataset = simu.data(n, p, p0, rho, nrep, SNR, seed)
x = dataset$x
y = dataset$y
beta = dataset$beta
mu = x%*%beta
sigma = dataset$sigma

## ------------------------------------------------------------------------
print(beta)

## ------------------------------------------------------------------------
library(BOSSreg)
# Choose a single replication as illustration
rep = seed
# Fit the model
boss_model = boss(x, y[,rep], intercept = FALSE)

## ------------------------------------------------------------------------
betahat_boss = boss_model$beta_boss
betahat_fs = boss_model$beta_fs
print(dim(betahat_boss))

## ---- fig.width=3, fig.height=3, fig.show='hold'-------------------------
# The heuristic degrees of freedom
plot(0:p, boss_model$hdf, main='hdf', ylab='', xlab='subset size', type='b')
abline(0, 1, lty=2)
# AICc-hdf (scaled by 1/n, and up to a constant)
plot(0:p, boss_model$IC_boss$aicc, main='AICc-hdf', ylab='', xlab='subset size', type='b')

## ------------------------------------------------------------------------
# The default is chosen by AICc
betahat_aicc = coef(boss_model)
muhat_aicc = predict(boss_model, newx=x)
# Use Cp rather than AICc
betahat_cp = coef(boss_model, ic='cp')
muhat_cp = predict(boss_model, newx=x)

## ------------------------------------------------------------------------
# The default is 10-fold CV with 1 replication
set.seed(seed)
boss_cv_model = cv.boss(x, y[,rep], intercept=FALSE)
# Coefficient vector selected by minimizing CV error
betahat_cv = coef(boss_cv_model)
# Fitted values
muhat_cv = predict(boss_cv_model, newx=x)

## ------------------------------------------------------------------------
# Coefficient vector for FS selected by CV
betahat_fs_cv = coef(boss_cv_model, method='fs')
# Fitted values
muhat_fs_cv = predict(boss_cv_model, newx=x, method='fs')

## ------------------------------------------------------------------------
tmp = cbind(betahat_aicc, betahat_cp, betahat_cv, betahat_fs_cv)
dimnames(tmp) = list(dimnames(tmp)[[1]], c('B0SS AICc', 'BOSS Cp', 'BOSS CV', 'FS CV'))
print(tmp)

## ------------------------------------------------------------------------
# X9 joins first
print(boss_model$steps_x)

## ---- eval=FALSE---------------------------------------------------------
#  # Function to calculate RMSE
#  calc.rmse <- function(muhat){
#    sqrt( Matrix::colSums(sweep(muhat, 1, mu)^2) / n )
#  }
#  rmse_solutionpath = list(BOSS=list(), FS=list())
#  for(rep in 1:nrep){
#    boss_model = boss(x, y[,rep], intercept=FALSE)
#    # RMSE along the solution path
#    rmse_solutionpath[['BOSS']][[rep]] = calc.rmse(x %*% boss_model$beta_boss)
#    rmse_solutionpath[['FS']][[rep]] = calc.rmse(x %*% boss_model$beta_fs)
#  }
#  # saveRDS(rmse_solutionpath, '/vignettes/rmse_solutionpath.rds')

## ---- include=FALSE------------------------------------------------------
rmse_solutionpath = readRDS(gzcon(url('https://raw.githubusercontent.com/sentian/BOSSreg/master/r-package/vignettes/rmse_solutionpath.rds')))

## ---- fig.width=5, fig.height=4------------------------------------------
# Average RMSE over replications
rmse_avg = lapply(rmse_solutionpath, function(xx){colMeans(do.call(rbind, xx))})
plot(0:p, rmse_avg$FS, col='blue', pch=1, main='Average RMSE along the solution path', 
     ylab='RMSE', xlab='Subset size')
points(0:p, rmse_avg$BOSS, col='red', pch=3)
legend('topright', legend = c('FS', 'BOSS'), col=c('blue', 'red'), pch=c(1,3))

## ---- eval=FALSE---------------------------------------------------------
#  rmse = nvar = list(BOSS=c(), FS=c())
#  set.seed(seed)
#  for(rep in 1:nrep){
#    boss_cv_model = cv.boss(x, y[,rep], intercept=FALSE)
#    # RMSE for the optimal subset selected via a selection rule
#    rmse[['BOSS']][rep] = calc.rmse(predict(boss_cv_model$boss, newx = x)) # AICc
#    rmse[['FS']][rep] = calc.rmse(predict(boss_cv_model, newx = x, method = 'fs')) # CV
#    # Number of variables
#    nvar[['BOSS']][rep] = sum(coef(boss_cv_model$boss)!=0)
#    nvar[['FS']][rep] = sum(coef(boss_cv_model, method='fs')!=0)
#  }
#  # saveRDS(list(rmse=rmse, nvar=nvar), '/vignettes/boss_fs.rds')

## ---- include=FALSE------------------------------------------------------
tmp = readRDS(gzcon(url('https://raw.githubusercontent.com/sentian/BOSSreg/master/r-package/vignettes/boss_fs.rds')))
rmse = tmp$rmse
nvar = tmp$nvar

## ---- fig.width=3, fig.height=3, fig.show='hold'-------------------------
# Make the plots
boxplot(rmse, outline=FALSE, main='RMSE')
boxplot(nvar, outline=FALSE, main='Number of predictors')

## ---- include=FALSE------------------------------------------------------
library(ISLR)
dataset = list()
# Boston Housing data
tmp = Boston
tmp = na.omit(tmp)
tmp$chas = as.factor(tmp$chas)
dataset$boston$x = data.matrix(tmp[,!names(tmp) %in% 'medv'])
dataset$boston$y = tmp$medv

# MLB hitters salary
tmp = Hitters
tmp = na.omit(tmp)
tmp[,c('League', 'Division', 'NewLeague')] = 
  lapply(tmp[,c('League', 'Division', 'NewLeague')], as.factor)
dataset$hitters$x = data.matrix(tmp[,!(names(tmp) %in% c('Salary'))])
dataset$hitters$y = tmp$Salary

# College data
tmp = College
tmp$Private = as.factor(tmp$Private)
dataset$college$x = data.matrix(tmp[,!(names(tmp) %in% c('Outstate'))])
dataset$college$y = tmp$Outstate

# Auto data
tmp = Auto
dataset$auto$x = data.matrix(tmp[,!(names(tmp) %in% c('mpg','name','origin'))])
dataset$auto$y = tmp$mpg

## ---- echo=FALSE---------------------------------------------------------
library(knitr)
library(kableExtra)
# read the results
result = readRDS(gzcon(url('https://raw.githubusercontent.com/sentian/BOSSreg/master/r-package/vignettes/realdata.rds')))
# function to extract the results
tmp_function <- function(method){
  unlist(lapply(result, function(xx){
    unlist(lapply(xx, function(yy){
      round(mean(yy[[method]]), 3)
      }))
    }))
}
tmp = data.frame(Dataset = rep(names(result), each=3),
  n_p = rep(unlist(lapply(dataset, function(xx){paste(dim(xx$x), collapse = ', ')})) , each=3),
  Metrics = rep(c('RMSE', '# predictors', 'running time (s)'), length(result)),
  BOSS = tmp_function('boss'),
  FS = tmp_function('fs'),
  LASSO = tmp_function('lasso'),
  SparseNet = tmp_function('sparsenet'))
rownames(tmp) = NULL
colnames(tmp)[2] = 'n, p'
kable(tmp, align = "c") %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "middle")

## ---- eval=FALSE---------------------------------------------------------
#  library(ISLR)
#  dataset = list()
#  # Boston Housing data
#  tmp = Boston
#  tmp = na.omit(tmp)
#  tmp$chas = as.factor(tmp$chas)
#  dataset$boston$x = data.matrix(tmp[,!names(tmp) %in% 'medv'])
#  dataset$boston$y = tmp$medv
#  
#  # MLB hitters salary
#  tmp = Hitters
#  tmp = na.omit(tmp)
#  tmp[,c('League', 'Division', 'NewLeague')] =
#    lapply(tmp[,c('League', 'Division', 'NewLeague')], as.factor)
#  dataset$hitters$x = data.matrix(tmp[,!(names(tmp) %in% c('Salary'))])
#  dataset$hitters$y = tmp$Salary
#  
#  # College data
#  tmp = College
#  tmp$Private = as.factor(tmp$Private)
#  dataset$college$x = data.matrix(tmp[,!(names(tmp) %in% c('Outstate'))])
#  dataset$college$y = tmp$Outstate
#  
#  # Auto data
#  tmp = Auto
#  dataset$auto$x = data.matrix(tmp[,!(names(tmp) %in% c('mpg','name','origin'))])
#  dataset$auto$y = tmp$mpg

## ---- eval=FALSE---------------------------------------------------------
#  library(glmnet)
#  library(sparsenet)
#  rmse <- function(y_hat, y){
#    sqrt(sum( (y_hat - y)^2 / length(y)) )
#  }
#  rdresult <- function(x, y, nrep, seed){
#    p = dim(x)[2]
#  
#    allmethods = c('lasso','sparsenet','boss','fs')
#    error = numvar = time = replicate(length(allmethods), rep(NA,nrep), simplify=F)
#    names(error) = names(numvar) = names(time) = allmethods
#  
#    set.seed(seed)
#    for(i in 1:nrep){
#      index = 1:nrow(x)
#      index = index[-i]
#  
#      x.train = x[index, , drop=FALSE]
#      y.train = y[index]
#      x.test = x[-index, , drop=FALSE]
#      x.test.withint = cbind(rep(1,nrow(x.test)), x.test)
#      y.test = y[-index]
#  
#      # BOSS
#      ptm = proc.time()
#      boss_model = boss(x.train, y.train, intercept = TRUE)
#      time_tmp = proc.time() - ptm
#      boss_pred = as.numeric( predict(boss_model, newx=x.test) )
#      error$boss[i] = rmse(boss_pred, y.test)
#      numvar$boss[i] = sum(coef(boss_model)!=0)
#      time$boss[i] = time_tmp[3]
#  
#      # FS
#      ptm = proc.time()
#      boss_cv_model = cv.boss(x.train, y.train)
#      time_tmp = proc.time() - ptm
#      fs_pred = as.numeric( predict(boss_cv_model, newx=x.test, method='fs') )
#      error$fs[i] = rmse(fs_pred, y.test)
#      numvar$fs[i] = sum(coef(boss_cv_model, method='fs')!=0)
#      time$fs[i] = time_tmp[3]
#  
#      # LASSO
#      ptm = proc.time()
#      lasso_model = glmnet(x.train, y.train, intercept=TRUE)
#      lasso_aicc = as.numeric(calc.ic(predict(lasso_model, newx=x.train), y.train,
#                                      ic='aicc', df=lasso_model$df+1))
#      lasso_pred = predict(lasso_model, newx=x.test, s=lasso_model$lambda[which.min(lasso_aicc)])
#      time_tmp = proc.time() - ptm
#      error$lasso[i] = rmse(lasso_pred, y.test)
#      numvar$lasso[i] = sum(coef(lasso_model, s=lasso_model$lambda[which.min(lasso_aicc)])!=0)
#      time$lasso[i] = time_tmp[3]
#  
#      # SparseNet
#      ptm = proc.time()
#      sparsenet_cv_model = cv.sparsenet(x.train, y.train)
#      time_tmp = proc.time() - ptm
#      sparsenet_pred = predict(sparsenet_cv_model, newx=x.test, which='parms.min')
#      error$sparsenet[i] = rmse(sparsenet_pred, y.test)
#      numvar$sparsenet[i] = sum(coef(sparsenet_cv_model, which='parms.min')!=0)
#      time$sparsenet[i] = time_tmp[3]
#    }
#    return(list(error=error, numvar=numvar, time=time))
#  }
#  result = lapply(dataset, function(xx){rdresult(xx$x, xx$y, nrow(xx$x), seed)})
#  # saveRDS(result, '/vignettes/realdata.rds')

## ---- eval=FALSE---------------------------------------------------------
#  library(knitr)
#  library(kableExtra)
#  # Function to extract the results
#  tmp_function <- function(method){
#    unlist(lapply(result, function(xx){
#      unlist(lapply(xx, function(yy){
#        round(mean(yy[[method]]), 3)
#        }))
#      }))
#  }
#  tmp = data.frame(Dataset = rep(names(result), each=3),
#    n_p = rep(unlist(lapply(dataset, function(xx){paste(dim(xx$x), collapse = ', ')})) , each=3),
#    Metrics = rep(c('RMSE', '# predictors', 'running time (s)'), length(result)),
#    BOSS = tmp_function('boss'),
#    FS = tmp_function('fs'),
#    LASSO = tmp_function('lasso'),
#    SparseNet = tmp_function('sparsenet'))
#  rownames(tmp) = NULL
#  colnames(tmp)[2] = 'n, p'
#  kable(tmp, align = "c") %>%
#    kable_styling(full_width = F) %>%
#    column_spec(1, bold = T) %>%
#    collapse_rows(columns = 1:2, valign = "middle")

