globalVariables(c("value", "variable", "terms", "strata"))
#' Boosting core function
#'
#' This function allows you to use gradient boosting for variable selection.
#' @param formula a formula object with a response value using the Surv
#'   function.
#' @param data a data.frame containing all variables specified in the formula.
#' @param rate the desired update rate used in the boosting algorithm.
#' @param num_iter an integer used as the number of iterations of the boosting
#'   algorithm. Default value is 500.
#' @param control_method specifies stopping method, options include: cv,
#'   num_selected, likelihood, BIC, AIC. Default is NULL, which will use a fixed
#'   number of iterations as specified by num_iter.
#' @param control_parameter is a list with the parameter(s) needed for each 
#' corresponding control_method option, the options are "cv_folds", "early_stop", 
#' "EBIC_gamma", "num_select", and "likelihood_tol." For cv method "cv_folds" 
#' specifies the number of cross validation folds (default is 10). For EBIC and 
#' AIC methods, "early_stop" is a TRUE/FALSE value for early stopping (default is FALSE). 
#' An additional parameter for the EBIC method is "EBIC_gamma" that is used to 
#' specify the penalty term, should be a value between 0 and 1. If using 
#' num_selected method, "num_select" will be the desired number of variables to 
#' select, should be an integer. If using likelihood as the method, 
#' "likelihood_tol" will be the small change in likelihood in which to stop 
#' once reached (default is 0.001). 
#' @param censoring_type currently only right censoring is implemented. 
#' @return a list containing the vector of coefficients ("beta"), variable
#'   selection matrix that contains the coefficients at each iteration
#'   ("selection_df"), the number of boosting iterations ("mstop"), and other
#'   stopping criteria if applicable to selected method. If using method BIC or AIC, the information criteria for each iteration is returned as a vector ("Information Criteria"). If using cross validation for stopping the criteria used for stopping is returned as a numeric vector ("cvrisk"). 
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' formula <- as.formula("Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + 
#' V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10" )
#' boosting_core(formula, data, rate=0.1, num_iter=500)
#' boosting_core(formula, data, rate=0.1, control_method="num_selected",
#' control_parameter=list(num_select = 5))
#' 
boosting_core <- function(formula, data, rate, num_iter=500, control_method=NULL, control_parameter=NULL, censoring_type = "right" ){
  # require(survival)  
  df <- stats::get_all_vars(formula, data=data) # takes a few seconds
  output <- list(formula=formula, data=data, rate=rate, censoring_type = censoring_type)
  data_name <- deparse(substitute(data))
  
  Call <- match.call()
  
  special <- c("strata")
  indx <- match(c("formula", "data", "na.action"), names(Call), nomatch = 0)
  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  
  temp$formula <- stats::terms(formula, special, data = data) # slightly slow
  coxenv <- new.env(parent = environment(formula))
  environment(temp$formula) <- coxenv
  mf <- eval(temp, parent.frame())
  Terms <- stats::terms(mf)
  
  # edits for Y
  Y <- stats::model.extract(mf, "response")
  time <- Y[,1]
  delta <- Y[,2]
  
  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- formula[1:2]
    xtemp <- formula[-2]
    #if (any(!is.na(match(xtemp, ytemp))))
    #  warning("a variable appears on both the left and right sides of the formula")
  }
  contrast.arg <- NULL
  stemp<- survival::untangle.specials(Terms, "strata", 1)
  if (length(stemp$vars) > 0) {
    dropterms <- stemp$terms
    temppred <- attr(terms, "predvars")
    Terms2 <- Terms[-dropterms]
    if (!is.null(temppred)) {
      attr(Terms2, "predvars") <- temppred[-(1 + dropterms)]
    }
    attr(Terms2,"intercept") <- 0
    X <- stats::model.matrix(Terms2, mf, constrasts = contrast.arg)
    renumber <- match(colnames(attr(Terms2, "factors")),
                      colnames(attr(Terms, "factors")))
    attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
  }
  else  X <- stats::model.matrix(Terms, mf, contrasts = contrast.arg)
  
  strats <- attr(Terms, "specials")$strata
  if (length(strats)) {
    stemp <- survival::untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) == 1)
      strata.keep <- mf[[stemp$vars]]
    else strata.keep <- survival::strata(mf[, stemp$vars], shortlabel = TRUE)
    strats <- as.numeric(strata.keep)
  }
  if(is.null(strats)){
    strats <- rep(1, length(delta))
    X_cov <- as.matrix(df[,-c(1,2)])
  }
  else{
    X_cov <- as.matrix(df[,-c(1,2,3)]) # need to also remove strata from X if it exists in df
  }
  
  adj_variables <- 0
  sample <- c(0:(length(delta)-1))
  num_strata = length(unique(strats))
  
  output <- list(data=data, rate=rate, censoring_type = censoring_type)
  # specify default values for control_parameter corresponding to each control_method
  if(is.null(control_parameter)){
    control_parameter = list(cv_folds = 10,
                             early_stop = FALSE,
                             EBIC_gamma = 0,
                             num_select = 1L,
                             likelihood_tol = 0.001)
  }
    
  # if want to input m_stop
  if(!is.null(num_iter) & is.null(control_method)){
    result <- boosting_stratify_core(sample, delta, strats, num_strata, X_cov, num_iter, rate, adj_variables)
    output$selection_df <- result
    output$beta <- result[num_iter,]
    output$mstop = num_iter
    control_method = "input"
  }
  # if cross validation
  else if(control_method == "cv"){
    if(is.null(control_parameter$cv_folds)){
      control_parameter$cv_folds=10
    }
    cv_result <- cross_validation_func_update(control_parameter$cv_folds, time, delta, X_cov, strats, rate, M_stop=num_iter)
    mstop_cv <- cv_result$mstop
    result <- boosting_stratify_core(sample, delta, strats, num_strata, X_cov, mstop_cv, rate, adj_variables)
    output$selection_df <- result
    output$beta <- result[mstop_cv,]
    output$mstop = mstop_cv
    output$cvrisk <- cv_result$cvrisk
    
  }
  # if # selected specified
  else if(control_method == "num_selected"){
    if(is.null(control_parameter$num_select)){
      control_parameter$num_select = 1
    }
    temp = boosting_stratify_numselected1(sample,  delta, strats, num_strata, X_cov, num_selected=control_parameter$num_select, rate, adj_variables)
    output$beta <- temp$beta
    output$mstop = temp$num_iterations
    output$selection_df <- temp$selection_df[1:temp$num_iterations,]
  }
  # if change in likelihood is specified
  else if(control_method == "likelihood"){
    if(is.null(control_parameter$likelihood_tol)){
      control_parameter$likelihood_tol = 0.001
    }
    temp = boosting_stratify_likelihood1(sample, delta, strats, num_strata, X_cov, rate, delta_likelihood=control_parameter$likelihood_tol, adj_variables)
    output$beta <- temp$beta
    output$mstop <- temp$num_iterations
    output$selection_df <- temp$selection_df
  }
  # if BIC is specified
  else if(control_method == "BIC"){
    if(is.null(control_parameter$early_stop)){
      control_parameter$early_stop = FALSE
    }
    if(is.null(control_parameter$EBIC_gamma)){
      control_parameter$EBIC_gamma = 0
    }
    temp = boosting_stratify_BIC1(sample, delta, strats, num_strata, X_cov, rate, early_stop=control_parameter$early_stop, adj_variables, gamma = control_parameter$EBIC_gamma)
    output$beta <- temp$beta
    output$mstop <- temp$num_iterations
    output$selection_df <- temp$selection_df[1:temp$num_iterations,]
    output$information.criteria <- temp$`Information Criteria`
  }
  else if(control_method == "AIC"){
    if(is.null(control_parameter$early_stop)){
      control_parameter$early_stop = FALSE
    }
    temp = boosting_stratify_BIC1(sample, delta, strats, num_strata, X_cov, rate, early_stop=control_parameter$early_stop, adj_variables, gamma=0, aic=TRUE)
    output$beta <- temp$beta
    output$mstop <- temp$num_iterations
    output$selection_df <- temp$selection_df[1:temp$num_iterations,]
    output$information.criteria <- temp$`Information Criteria`
  }
  else{
    cat("Incorrect option for mstop_method. Options are: 'cv', 'num_selected', 'likelihood', 'AIC', 'BIC', or provide a number for num_iter." )
    # break
  }
  output$control_method = control_method
  output$coefficients <- as.vector(output$beta)
  output$strata <- strats
  names(output$strata) <- names(strats)
  if("(Intercept)" %in% colnames(X)){
    names(output$coefficients) <- colnames(X)[2:ncol(X)]
  }
  else{
    names(output$coefficients) <- colnames(X)
  }
  output$formula <- formula
  output$call <- Call
  output$delta <- delta
  names(output)[[1]] <- data_name
  class(output) <- "boosting"
  output
}

#' Prints the call and coefficients from boosting model selection
#'
#' This function displays the coefficient estimates of all variables from a
#' model generated with the boosting_core function.
#' @param x output from boosting_core function.
#' @param ... ignored
#' @return list containing the coefficient vector and function call.
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' formula <- as.formula("Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + 
#' V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10" )
#' boosting.output <- boosting_core(formula, data, rate=0.1, num_iter=500)
#' print(boosting.output)
#' 
print.boosting = function(x, ...){
  x = x[c("call", "coefficients")]
  NextMethod()
}

#' Summary of boosting model selection
#'
#' This function displays the variables selected from a model generated with the
#' boosting_core function.
#' @param object output from boosting_core function.
#' @param all_beta default value is FALSE. If this is set to TRUE the 
#' coefficient estimates for all the parameters will be printed. 
#' @param ... ignored
#' @return list containing the coefficient vector, number of boosting
#'   iterations, and resulting formula from the variable selection.
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' formula <- as.formula("Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + 
#' V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10" )
#' boosting.output <- boosting_core(formula, data, rate=0.1, num_iter=500)
#' summary.boosting(boosting.output)
#' 
summary.boosting <- function(object, all_beta = NULL, ...)
{
  if(!is.null(cl <- object$call)){
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  
  cat("\n", "data: ", names(object)[1], "\n") 
  cat("\n", "n = ", dim(object[[1]])[1],"\n", "Number of events = ", sum(object$delta),"\n",
      #   "Number of strata = ", length(unique(x$data["strata"])), "\n",
      "Number of boosting iterations: mstop = ", object$mstop, "\n", "Step size = ", object$rate)
  cat("\n")
  if(is.null(all_beta)){
    reduced_beta <- object$coefficients[which(object$coefficients!=0)]
    cat("\n", "Coefficients:", "\n")
    print(reduced_beta)
    output <- list(reduced_beta, object$mstop)
  }
  else{
    cat("\n", "Coefficients:", "\n")
    print(object$coefficients)
  }

  # generate output but return as invisible 
  df_coxph <- object[[1]]
  df_coxph$strata <- object$strata
  strata.col <- which(sapply(object[[1]], identical, object$strata)==TRUE)
  strata.name <- names(object[[1]])[strata.col]
  strata.fmla <- paste("strata(",strata.name ,")")
  fmla_reduced <- stats::as.formula(paste("Surv(time,delta) ~ ", paste(c(names(reduced_beta), strata.fmla), collapse= "+")))
  # cat("\n", "Formula: ", deparse(fmla_reduced), "\n")
  output$formula <- fmla_reduced
  names(output) <- c("Coefficients", "Number of iterations", "Formula")
  invisible(output) # return(output) 
}

#' Boosting plot function
#'
#' This function allows you to visualize the coefficient paths of the boosting
#' algorithm.
#' @param x output from the boosting_core function.
#' @param y y coordinates of plot, default is NULL. 
#' @param type specifies type of coefficient plot. Default value is frequency
#' which plots the proportion of variables selected. Alternatively type set to 
#' "coefficients" plots the coefficient path for each variable.
#' @param ... ignored 
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' formula <- as.formula("Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + 
#' V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10" )
#' boosting.output <- boosting_core(formula, data, rate=0.1, num_iter=500)
#' plot.boosting(boosting.output)
#' plot.boosting(boosting.output, type="coefficients")
#' 
plot.boosting <- function(x, y=NULL, type="frequency", ...)
{
  selection_df <- x$selection_df
  type.coefficients <- pmatch("coef", type) 
  
  if(type == "frequency"){
    selection_df[which(selection_df!=0)] <- 1
    y <- rowMeans(selection_df)
    x_axis <- c(1:x$mstop)
    plot_data <- data.frame(cbind(x_axis,y))
    
    ggplot2::ggplot(plot_data, ggplot2::aes(x_axis)) +
      ggplot2::geom_line(ggplot2::aes(y = y, colour = "y"))  + ggplot2::theme_bw() + ggplot2::theme(legend.position="none", text=ggplot2::element_text(size=16)) +
      ggplot2::xlab("Number of boosting iterations") + ggplot2::ylab("Proportion of variables selected")  
  } else if(!is.na(type.coefficients)){
    selection_df <- rbind(rep(0,ncol(selection_df)), selection_df)
    x_axis <- c(0:x$mstop)
    plot_data <- data.frame(x_axis, selection_df)
    long <- reshape2::melt(plot_data, id.vars = c("x_axis"))
    ggplot2::ggplot(long, ggplot2::aes(x=x_axis, y=value, group=variable)) +
      ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::theme( text= ggplot2::element_text(size=16)) +
      ggplot2::ylab("Coefficient Estimate") + ggplot2::xlab("Number of iterations") +
      directlabels::geom_dl(ggplot2::aes(label = variable), method = list(directlabels::dl.combine("last.points"), cex = 1.2))
  }
}

#' Boosting inference function
#'
#' This function provides post selection inference.
#' @param x output from boosting_core function. 
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' formula <- as.formula("Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + 
#' V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10" )
#' boosting.output <- boosting_core(formula, data, rate=0.1, num_iter=500)
#' inference.boosting(boosting.output)
#'
# POST SELECTION INFERENCE
inference.boosting <- function(x){ 
  data = data.frame(x$data) 
  data$strata = x$strata
  variables.selected <- names(x$coefficients)[which(abs(x$coefficients) > x$rate)]
  if(!is.null(x$strata)){
    fmla = stats::formula(paste(x$formula[2] , "~ ",   paste(c(variables.selected, "strata(strata)"), collapse = "+")))
  }else{
    fmla = stats::formula(paste(x$formula[2] , "~ ", paste(variables.selected, collapse = " + ")))
  }
  fit <- survival::coxph(fmla, data=data)
  print(summary(fit))
}

#' Boosting predict function
#'
#' This function predicts the hazard ratio for each subject in the input dataset. 
#' @param object output from boosting_core function.
#' @param newdata data.frame used for prediction. Default is NULL and will 
#' use data specified for boosting algorithm.
#' @param ... ignored
#' @return vector of the hazard ratio for each observation relative to the 
#' sample average.
#' @keywords gradient boosting
#' @export
#' @examples
#' data <- simulate_survival_cox(true_beta=c(1,1,1,1,1,0,0,0,0,0))
#' formula <- as.formula("Surv(time,delta) ~ strata(strata_idx) + V1 + V2 + 
#' V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10" )
#' boosting.output <- boosting_core(formula, data, rate=0.1, num_iter=500)
#' predict.boosting(boosting.output)
#'
# PREDICT BOOSTING
predict.boosting <- function(object, newdata=NULL, ...){
  if(is.null(newdata)) newdata <- object$data
  
  var_selected <- names(which(object$coefficients!=0))
  sample_avg <- colMeans(newdata[, var_selected])
  
  #calculate exp(Xb) for sample average
  hazard_avg <- exp(sample_avg %*% object$coefficients[which(object$coefficients!=0)])
  
  predict <- NULL
  for(i in 1:nrow(newdata)){
    #calculate exp(X'b) for individual
    current_values <- unlist(newdata[i,var_selected], use.names = FALSE)
    hazard_ind <- exp(current_values %*% object$coefficients[which(object$coefficients!=0)])
    predict[i] <- hazard_ind/hazard_avg
  }
  return(predict)
}