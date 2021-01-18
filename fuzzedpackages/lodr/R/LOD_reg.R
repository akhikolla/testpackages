#################                     Function
lod_lm <- function(data, frmla, lod=NULL, var_LOD=NULL,
                           nSamples=250,
                           fill_in_method="mean",
                           convergenceCriterion=0.001,
                           boots=25){
  cl <- match.call()
  finalEstimates <- list()
  
######################################################################
# Check if lod/var_LOD are null: Return lm() fn call                 #
######################################################################

if(is.null(lod)==TRUE|is.null(var_LOD)==TRUE){
  warning("No LOD covariates and/or no LOD values provided: returning lm() call output")
  return(lm(frmla, data))
}
  
if(all(var_LOD%in%colnames(data))==FALSE){
  stop("LOD covariate(s) not in dataset")
}
  
if(is.numeric(lod)==FALSE|any(is.na(lod))==TRUE){
  stop("Argument lod not specified correctly; must be numeric vector with no missing values")
}
  
if(is.numeric(nSamples)==FALSE|length(nSamples)>1|all(nSamples==as.integer(nSamples))==FALSE){
  stop("Argument nSamples must be integer")
}
  
if(is.numeric(boots)==FALSE|length(boots)>1|all(boots==as.integer(boots))==FALSE){
  stop("Argument boots must be integer")
}
  
if(is.numeric(convergenceCriterion)==FALSE|length(convergenceCriterion)>1|convergenceCriterion<0){
  stop("Argument convergenceCriterion must be positive number")
}
  
######################################################################
# Check if boots>0, else print "Boots=0, SEs not computed            #
######################################################################

if(boots<0){
  stop("boots must be a nonnegative integer")
}
  
if(boots==0){
  message("NOTE: boots=0 specified, SEs not computed and marked by NA")
}
  
######################################################################
# Create required datasets                                           #
######################################################################

dataset_Y <- model.frame(frmla,data=data)[1]
dataset_X <- data.frame(model.matrix(frmla, data))[,-1]

# Reorder dataset_X so LOD vars are last
dataset_X_reorder <- 
  data.frame(dataset_X[,!names(dataset_X)%in%var_LOD],
             dataset_X[,names(dataset_X)%in%var_LOD])
names(dataset_X_reorder)[!names(dataset_X)%in%var_LOD] <- 
  names(dataset_X)[!names(dataset_X)%in%var_LOD]

names(dataset_X_reorder)[names(dataset_X)%in%var_LOD] <- 
  names(dataset_X)[names(dataset_X)%in%var_LOD]

subData <- data.frame(cbind(dataset_Y,Intercept=1,dataset_X_reorder))

Data <- subData
sub2Data <- subData

sqrt2LOd <- function(x){
  ifelse(x>0, x/sqrt(2), x*sqrt(2))
}
subd <- sqrt2LOd(lod)

for(i in 1:length(var_LOD)){
  Data[[var_LOD[i]]] <- ifelse(subData[[var_LOD[i]]]>lod[i],
                               subData[[var_LOD[i]]], NA)
  
  sub2Data[[var_LOD[i]]] <- ifelse(subData[[var_LOD[i]]]>lod[i],
                                   subData[[var_LOD[i]]], subd[i])
}

# Create data for analysis
## Obs (Data), complete case (ccData), sub (subData), and sub sqrt2 (sub2Data) datasets
ccData <- Data[complete.cases(Data),]

######################################################################
# Enter convergence criterion, LOD, etc.                             #
######################################################################

n <- dim(Data)[1]
nObservations <- n

#######################################################################
# Perform Complete-Case, Substitution Analysis                        #
#######################################################################

ccModel <- glm( frmla,
                family = gaussian(), data = ccData )

######################################################################
# Perform all Substitution Analyses                                  #
######################################################################

# Using LOD for sub
ccSub_LOD <- glm( frmla,
                  family = gaussian(), data = subData )

# Using LOD/sqrt(2) for sub
ccSub_LODsqrt2 <- glm( frmla,
                       family = gaussian(), data = sub2Data )

######################################################################
# Obtain parameter estimates to be used as initial estimates in ARMS #
######################################################################

# Extract Beta and residual variance
BetaEstimatesCC <- as.numeric( summary( ccModel )$coefficients[,1])
names(BetaEstimatesCC) <- rownames(summary( ccModel )$coefficients)
BetaEstimatesSubsqrt2 <- as.numeric( summary( ccSub_LODsqrt2 )$coefficients[,1])
names(BetaEstimatesSubsqrt2) <- rownames(summary( ccSub_LODsqrt2 )$coefficients)

# Reorder Beta so LOD vars are last
beta_names_reordered <- 
  c(names(BetaEstimatesCC)[!names(BetaEstimatesCC)%in%var_LOD],
    names(BetaEstimatesCC)[names(BetaEstimatesCC)%in%var_LOD])
BetaEstimatesCC_reorder <- BetaEstimatesCC[beta_names_reordered]
BetaEstimatesSubsqrt2_reorder <- BetaEstimatesSubsqrt2[beta_names_reordered]

# Extract mean vector estimate and covariance matrix estimate of covariates
cat_var <- names(ccModel$contrasts)
length_unique <- function(x){length(unique(x))}
unique_values <- apply(dataset_X_reorder, MARGIN=2,FUN=length_unique)
binary_vars <- names(unique_values[unique_values<3])
remove_vars <- c(var_LOD, cat_var, binary_vars)
var_noLOD <- names(dataset_X_reorder[,!(names(dataset_X_reorder)%in%remove_vars),
                                     drop=FALSE])
var_keep <- names(ccData)[names(ccData) %in% c(var_noLOD, var_LOD)]

# Reorder dataset before calcs so LOD vars are last
xMeanCC <- apply( ccData[,var_keep], 2, mean ) #-1 to eliminate outcome Y
xCovCC <- cov( ccData[,var_keep] )

#######################################################################
# Perform ARMS MLE Sampling using C++ compiled code                                          #
#######################################################################

## Create matrix to hold limits of detection

LOD_mat <- cbind(rep(-100, dim(Data[,-1])[2]), rep(NA, dim(Data[,-1])[2]))
LOD_mat[which(names(Data[,-1])%in%var_LOD),2] <- lod

# Estimation
est_obj <- LOD_fit(y_data=Data[,1], 
                   x_data=as.matrix(Data[,-1]),
                   mean_x_preds=xMeanCC,
                   beta=BetaEstimatesCC_reorder,
                   sigma_2_y = sigma(ccModel)^2,
                   sigma_x_preds = xCovCC,
                   no_of_samples=nSamples, 
                   threshold = convergenceCriterion, 
                   max_iterations = 100,
                   LOD_u_l = LOD_mat)

LOD_ests <- est_obj$beta_estimate_last_iteration
names(LOD_ests) <- names(BetaEstimatesCC_reorder)

# Bootstrap SEs
if(boots>0){
boot_obj <- LOD_bootstrap_fit(num_of_boots=boots,
                              y_data=Data[,1], 
                              x_data=as.matrix(Data[,-1]),
                              no_of_samples=nSamples, 
                              threshold = convergenceCriterion, 
                              max_iterations = 100,
                              LOD_u_l = LOD_mat)
boot_SE_reorder <- apply(do.call("rbind", boot_obj), 2, sd)

}else{
  boot_SE_reorder <- rep(NA, length(BetaEstimatesCC_reorder))
}
names(boot_SE_reorder) <- names(BetaEstimatesCC_reorder)

# Create lm_lod object
finalEstimates$coefficients <- 
  LOD_ests[colnames(model.matrix(frmla, data))]
finalEstimates$boot_SE <- 
  boot_SE_reorder[colnames(model.matrix(frmla, data))]
finalEstimates$rank <- dim(Data[,-1])[2]
finalEstimates$residuals <- Data[,1]-finalEstimates$fitted.values
finalEstimates$df.residual <- n-dim(Data[,-1])[2]
if(!is.null(cat_var)){
  finalEstimates$xlevels <- list()
  for(i in 1:length(cat_var)){
    finalEstimates$xlevels[[cat_var]] <- levels(data[[cat_var]])
  }
}

# Fill in values outside of LOD using method of choice, then calc fitted values and residuals
fill_in_data <- Data
if(fill_in_method=="mean"){
  for(i in 1:length(var_LOD)){
    fill_in_data[is.na(fill_in_data[,var_LOD[i]]),var_LOD[i]] <- 
      mean(data[,var_LOD[i]], na.rm = TRUE)
  }
}else{
  if(fill_in_method=="LOD"){
  for(i in 1:length(var_LOD)){
    fill_in_data[is.na(fill_in_data[,var_LOD[i]]),var_LOD[i]] <- lod[i]
  }
  }else{
      stop("Argument method must be = mean or = LOD")
    }
}

finalEstimates$fitted.values <- as.matrix(fill_in_data[,-1])%*%LOD_ests
finalEstimates$residuals <- fill_in_data[,1]-finalEstimates$fitted.values
finalEstimates$model <- Data[,names(model.frame(frmla,data=data))]
finalEstimates$terms <- 
  glm( frmla,
       family = gaussian(), data = data )$terms
finalEstimates$call <- cl

class(finalEstimates) <- "lod_lm"
return(finalEstimates)
}
##############################################################

## Create new generic fns: summary, print, coef, effects, residuals, fitted, vcov
# print
print.lod_lm <- function(x,...){
  print(list("call"=x$call,
             "coefficients"=x$coefficients))
}

# summary
summary.lod_lm <- function(object,...){
    output_obj <- list()
    param_values <- as.list(object$call)
  
    # coefficients
    coefficients_mat <- matrix(nrow=dim(model.matrix(eval(param_values$frmla),
                                                    eval(param_values$data)))[2],
                               ncol=4)
    rownames(coefficients_mat) <- colnames(model.matrix(eval(param_values$frmla),
                                                        eval(param_values$data)))
    colnames(coefficients_mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    coefficients_mat[,"Estimate"] <- object$coefficients
    coefficients_mat[,"Std. Error"] <- object$boot_SE
    coefficients_mat[,"t value"] <- coefficients_mat[,"Estimate"]/coefficients_mat[,"Std. Error"]
    coefficients_mat[,"Pr(>|t|)"] <- 2*(1-pt(abs(coefficients_mat[,"t value"]),
                                          df=dim(model.matrix(eval(param_values$frmla),
                                                              eval(param_values$data)))[1]-
                                            dim(model.matrix(eval(param_values$frmla),
                                                              eval(param_values$data)))[2]))
    output_obj$coefficients <- coefficients_mat
    
    # add in other components found in summary.lm, expect covariance matrix of coef estimates
    output_obj$call <- object$call
    output_obj$residuals <- object$residuals
    output_obj$df <- c(object$rank, object$df.residual)
    output_obj$sigma <- sum((object$residuals)^2)/(object$df.residual)
    # ssm <- sum((object$fitted.values-mean(object$model[,1]))^2)
    # f_stat <- 
    #   (ssm/(object$rank-1))/output_obj$sigma
    # output_obj$fstatistic <- c(f_stat, object$rank-1, object$df.residual)
    # output_obj$r.squared <- 
    #   1-(sum((object$residuals)^2)/sum((object$model[,1]-mean(object$model[,1]))^2))
    # output_obj$adj.r.squared <-
    #   1-(1-output_obj$r.squared)*((dim(object$mode)[1]-1)/(object$df.residual-1))
    class(output_obj) <- "summary.lod_lm"
    
    return(output_obj)
  }
  
# print summary
print.summary.lod_lm <- function(x,...){
  print(x$coefficients)
}

# coef
coef.lod_lm <- function(object,...){
  print(object$coefficients)
}

# residuals
residuals.lod_lm <- function(object,...){
  print(object$residuals)
}

# fitted.values
fitted.lod_lm <- function(object,...){
  print(object$fitted.values)
}
