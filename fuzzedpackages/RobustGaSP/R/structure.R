##########################################################################
## Class definition
## 
## Robust GaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2015-present Mengyang Gu, Jesus Palomo , James O. Berger
##							  
##    
##########################################################################

## rgasp Class

setClass("rgasp", 		
         representation( 
           p = "integer",          ## dimension of the inputs
           num_obs = "integer",          ## observations number
           ## data
           input = "matrix",           ## the design of experiments, size nxp
           output = "matrix",           ## the observations, size nx1
           X="matrix",                  ## mean basis, size nxq
           zero_mean="character",      ## Yes or no
           q="integer",                 ## number of mean basis
           LB="vector",                  ## lower bound for inverse range parameters beta px1
           beta_initial="vector",       ###initial values of inverse range parameters px1
           beta_hat="vector",                ###estimated inverse range parameters px1
           log_post="numeric",           ####logarithm of marginal posterior
           R0="list",                    ##abs difference of each type of input
          ###R="matrix",                  ## correlation matrix, size nxn
           theta_hat="vector",          ## regression parameter for mean, size qx1
           L="matrix",                  ####cholesky decomposition correlation matrix R, nxn, L%*%t(L)=R
          # R.inv_X="matrix" ,           ####R.inv%*%X, to avoid compute it again for prediction
          # Xt_R.inv_X.inv="matrix",     ####(t(X)%*%R.inv%*%X)^{-1}, to avoid compute it again for prediction     
           sigma2_hat="numeric",                 ####hat sigma^2
           LX="matrix",                  ####cholesky decomposition of correlation matrix t(X)%*%R^{-1}%*%X, qxq
           CL="vector",                  ##### used for lower bound and the prior
           nugget="numeric",                  #####nugget variance ratio 
           nugget.est="logical",         ### nugget is estimated or fixed
           kernel_type="vector",       #####type of kernel to specify
           alpha="vector",                 ####roughness parameter in the kernel
           method="character",          ##### post_mode, MLE, MMLE
           isotropic="logical",          ##isotropic or separable kernel
           call = "language"         ## user call
         ), 
         validity = function(object) {
           if (object@num_obs <= object@p) {
             return("the number of experiments must be larger than the spatial dimension")
           }
           
           if (ncol(object@output) != 1) {
             return("the response must have one dimension")
           }
           
           if (!identical(nrow(object@input), nrow(object@output))) {
             return("the number of observations is not equal to the number of experiments")
           }
           
           TRUE
         }
)
if(!isGeneric("show")) {
  setGeneric(name = "show",
             def = function(object) standardGeneric("show")
  )
}

setMethod("show", "rgasp",
          function(object){
            show.rgasp(object)		
          }
)

setClass("predrgasp",
         representation(
           call = "language",         ## user call
           mean = "numeric",          ## the mean value of the prediction
           lower95 = "numeric",       ## the lower 95 bound of the prediction
           upper95 = "numeric",       ## the upper 95 bound of the prediction
           sd = "numeric"             ## the standard deviation of the prediction
         )
)

if(!isGeneric("as.S3prediction")) {
  setGeneric(name = "as.S3prediction",
             def = function(object, ...) standardGeneric("as.S3prediction")
  )
}

setMethod("as.S3prediction", "predrgasp",
          definition = function(object)
          {
            structure(list(mean = object@mean, lower95 = object@lower95, 
                           upper95 = object@upper95, sd = object@sd), class = "predrgasp")
          }
)

if(!isGeneric("as.S4prediction")) {
  setGeneric(name = "as.S4prediction",
             def = function(object, ...) standardGeneric("as.S4prediction")
  )
}

setMethod("as.S4prediction", "predrgasp",
          definition = function(object)
          {
            as.S4prediction.predict(object=object)
          }
)


if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object, ...) standardGeneric("predict")
  )
}

setMethod("predict", "rgasp",
          definition=function(object, testing_input, testing_trend=matrix(1,dim(testing_input)[1],1),
                              r0=NA, 
                              interval_data=T,
                              outasS3 = T, ...) {
            predict.rgasp(object = object, testing_input = testing_input, 
                          testing_trend=testing_trend , r0=r0, 
                          interval_data=T,outasS3 = outasS3,...)
          }
)




if(!isGeneric("simulate")) {
  setGeneric(name = "simulate",
             def = function(object, ...) standardGeneric("simulate")
  )
}

setMethod("simulate", "rgasp",
          definition=function(object, testing_input, num_sample=1,testing_trend=matrix(1,dim(testing_input)[1],1), 
                              r0=NA, rr0=NA,
                              sample_data=T,...) {
            simulate.rgasp(object = object, testing_input = testing_input, num_sample=num_sample,
                           testing_trend=testing_trend,
                           r0=r0, rr0=rr0,sample_data=sample_data, ...)
          }
)

if(!isGeneric("Sample")) {
  setGeneric(name = "Sample",
             def = function(object, ...) standardGeneric("Sample")
  )
}

setMethod("Sample", "rgasp",
          definition=function(object, testing_input, num_sample=1,
                              testing_trend=matrix(1,dim(testing_input)[1],1), 
                              r0=NA, rr0=NA,sample_data=T,...) {
            simulate.rgasp(object = object, testing_input = testing_input, num_sample=num_sample,
                           testing_trend=testing_trend,r0=r0, rr0=rr0,sample_data=sample_data,
                           ...)
          }
)

# if(!isGeneric("Sample.rgasp")) {
#   setGeneric(name = "Sample.rgasp",
#              def = function(object, ...) standardGeneric("Sample.rgasp")
#   )
# }
# 
# setMethod("Sample.rgasp", "rgasp",
#           definition=function(object, testing_input, num_sample=1,testing_trend=matrix(1,dim(testing_input)[1],1), ...) {
#             simulate.rgasp(object = object, testing_input = testing_input, num_sample=num_sample,
#                            testing_trend=testing_trend , ...)
#           }
# )



if(!isGeneric("plot")) {
  setGeneric(name = "plot",
             def = function(x, y, ...) standardGeneric("plot")
  )
}

setMethod("plot", 
          signature(x = "rgasp"),
          definition=function(x, y, ...) {
            if (!missing(y)) warning("Argument y is ignored (not used in this plot method)")
            plot.rgasp(x=x, ...)
          
          }
)


##ppgasp
setClass("ppgasp", 		
         representation( 
           p = "integer",          ## dimension of the inputs
           num_obs = "integer",          ## observations number
           k = "integer",          ## locations number
           ## data
           input = "matrix",           ## the design of experiments, size nxp
           output = "matrix",           ## the observations, size nxk
           X="matrix",                  ## mean basis, size nxq
           zero_mean="character",      ## Yes or no
           q="integer",                 ## number of mean basis
           LB="vector",                  ## lower bound for inverse range parameters beta px1
           beta_initial="vector",       ###initial values of inverse range parameters px1
           beta_hat="vector",                ###estimated inverse range parameters px1
           log_post="numeric",           ####logarithm of marginal posterior
           R0="list",                    ##abs difference of each type of input
           ###R="matrix",                  ## correlation matrix, size nxn
           theta_hat="matrix",          ## regression parameter for mean, size qx1
           L="matrix",                  ####cholesky decomposition correlation matrix R, nxn, L%*%t(L)=R
           # R.inv_X="matrix" ,           ####R.inv%*%X, to avoid compute it again for prediction
           # Xt_R.inv_X.inv="matrix",     ####(t(X)%*%R.inv%*%X)^{-1}, to avoid compute it again for prediction     
           sigma2_hat="vector",                 ####hat sigma^2
           LX="matrix",                  ####cholesky decomposition of correlation matrix t(X)%*%R^{-1}%*%X, qxq
           CL="vector",                  ##### used for lower bound and the prior
           nugget="numeric",                  #####nugget variance ratio 
           nugget.est="logical",         ### nugget is estimated or fixed
           kernel_type="vector",       #####type of kernel to specify
           alpha="vector",                 ####roughness parameter in the kernel
           method="character",          ##### post_mode, MLE, MMLE
           isotropic="logical",          ##isotropic or separable kernel
           call = "language"         ## user call
         ), 
         validity = function(object) {
           if (object@num_obs <= object@p) {
             return("the number of experiments must be larger than the spatial dimension")
           }
           
           if (ncol(object@output) != 1) {
             return("the response must have one dimension")
           }
           
           if (!identical(nrow(object@input), nrow(object@output))) {
             return("the number of observations is not equal to the number of experiments")
           }
           
           TRUE
         }
)


setClass("predppgasp",
         representation(
           call = "language",         ## user call
           mean = "matrix",          ## the mean value of the prediction
           lower95 = "matrix",       ## the lower 95 bound of the prediction
           upper95 = "matrix",       ## the upper 95 bound of the prediction
           sd = "matrix"             ## the standard deviation of the prediction
         )
)

setMethod("as.S3prediction", "predppgasp",
          definition = function(object)
          {
            structure(list(mean = object@mean, lower95 = object@lower95, 
                           upper95 = object@upper95, sd = object@sd), class = "predppgasp")
          }
)

setMethod("as.S4prediction", "predppgasp",
          definition = function(object)
          {
            as.S4prediction.predict_ppgasp(object=object)
          }
)


setMethod("show", "ppgasp",
          function(object){
            show.ppgasp(object)		
          }
)


setMethod("predict", "ppgasp",
          definition=function(object, testing_input, testing_trend=matrix(1,dim(testing_input)[1],1),
                              r0=NA, 
                              interval_data=T,
                              outasS3 = T,...) {
            predict.ppgasp(object = object, testing_input = testing_input,r0=r0, 
                           interval_data=interval_data,
                          testing_trend=testing_trend, outasS3 = outasS3,...)
          }
)
