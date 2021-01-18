
#' @useDynLib SSLR



#' @title Predictions of model_sslr_fitted class
#' @description Predicts from model. There are different types: class, prob, raw
#' class returns tibble with one column
#' prob returns tibble with probabilities class columns
#' raw returns factor or numeric values
#' @param object  model_sslr_fitted model built.
#' @param x A object that can be coerced as matrix.
#' Depending on how was the model built, \code{x} is interpreted as a matrix
#' with the distances between the unseen instances and the selected training instances,
#' or a matrix of instances.
#' @param ... This parameter is included for compatibility reasons.
#' @param type of predict in principal model: class, raw, prob, vote, max_prob, numeric
#' @return tibble or vector.
#' @method predict model_sslr_fitted
#' @importFrom stats predict
#' @importFrom magrittr %>%
#' @importFrom dplyr as_tibble
#' @export
predict.model_sslr_fitted <- function(object, x,type = NULL,...){

  if(is.null(type)){
    if(object$model$mode == "classification")
      type <- "class"
    else
      type = "numeric"
  }

  if (!(type %in% object$model$pred.params))
    stop("No this type prediction module defined for this model.", call. = FALSE)

  #Apply formula if exists
  if(!is.null(object$formula)){
    x <- get_x_y(object$formula,x)$x
  }
  else{
    x <- x[,object$col_names_fit,drop = FALSE]
  }

  result <- object$model %>% predict(x, type = type)

  if(is.factor(result) & type == "class")
    pred_factor_tibble(result)

  else if(is.vector(result) & type == "numeric")
    pred_numeric_tibble(result)


  else if(type == 'prob'){
    result <- as.data.frame(result)
    names(result) <- paste0(".pred_", names(result))
    as_tibble(result)
  }

  else
    return(result)


}

#' @title Predictions of unlabeled data
#' @description Predictions of unlabeled data (transductive)
#' raw returns factor or numeric values
#' @param object  model_sslr_fitted model built
#' @param type of predict in principal model: class, raw
#' @param ... other parameters to be passed
#' @export predictions.model_sslr_fitted
#' @export
predictions.model_sslr_fitted <- function(object,type = "class",...){


  result <- object$model %>% predictions()

  if(type == "class")
    result <- pred_factor_tibble(result)

  return(result)
}




#' @title FUNCTION TO GET REAL X AND Y WITH FORMULA AND DATA
#' @description FUNCTION TO GET REAL X AND Y WITH FORMULA AND DATA
#' @param form formula
#' @param data data values, matrix, dataframe..
#' @return x (matrix,dataframe...) and y(factor)
#' @importFrom plyr is.formula
#' @import stats
get_x_y <- function(form, data) {

  composition <- class(data)

  if ("data.frame" %in% composition)
    composition <- "data.frame"

  if (!(composition %in% c("data.frame", "matrix")))
    stop("`composition` should be either 'data.frame' or ",
         "'matrix'.",
         call. = FALSE)

  formula <- as.formula(form)

  if (!is.formula(formula))
    stop("Formula is not a formula", call. = FALSE)

  col_name.class <- as.character(formula[2]) #name of class

  if (col_name.class == "" | is.null(col_name.class))
    stop("Not exists objective var", call. = FALSE)

  #EXECUTE MODEL FRAME
  datos <- stats::model.frame(formula = formula, data = as.data.frame(data),
                              na.action = NULL)

  #GET Y
  y <- datos[, c(col_name.class)]

  #GET X
  x <- datos
  x[, col_name.class] <- NULL

  #IF NOT DATAFRAME , RETURN MATRIX IN X
  if (composition == "data.frame") {
    x <- as.data.frame(x)
  }
  else {
    x <- as.matrix(x)
  }

  list(x = x, y = y)
}

#' @title FUNCTION TO GET FUNCTION METHOD
#' @description FUNCTION TO GET FUNCTION METHOD SPECIFIC
#' @param met character
#' @return  method_train (function)
get_function <- function(met) {
  method <- toupper(as.character(met))

  #print(paste("METHOD IS", method))


  #CHECK IF IT EXISTS
  if (!(method %in% c("COBC", "DEMOCRATIC", "SELFTRAINING", "SETRED", "SNNRCE", "TRITRAINING"))) {
    stop("This method not exist")
  }

  model_train <- NULL
  method_train <- NULL

  #ASSIGN FUNC

  if (method == "SELFTRAINING") {

    method_train <- selfTraining
  }

  else if (method == "COBC") {
    method_train <- coBC
  }

  else if (method == "DEMOCRATIC") {
    method_train <- democratic
  }

  else if (method == "SETRED") {
    method_train <- setred
  }

  else if (method == "SNNRCE") {
    method_train <- snnrce
  }

  else if (method == "TRITRAINING") {
    method_train <- triTraining
  }

  method_train

}



#' @title FUNCTION TO GET FUNCTION METHOD
#' @description FUNCTION TO GET FUNCTION METHOD GENERIC
#' @param met character
#' @return  method_train (function)
get_function_generic <- function(met) {
  #GET METHOD TO APPLY
  method <- toupper(as.character(met))

  #print(paste("METHOD IS", method))

  #CHECK IF IT EXISTS
  if (!(method %in% c("COBCG", "DEMOCRATICG", "SELFTRAININGG", "SETREDG", "SNNRCEG", "TRITRAININGG",
                     "COBC", "DEMOCRATIC", "SELFTRAINING", "SETRED", "SNNRCE", "TRITRAINING"))) {
    stop("This method not exist")
  }

  model_train <- NULL
  method_train <- NULL

  #ASSIGN FUNC

  if (method %in% c("SELFTRAININGG", "SELFTRAINING")) {

    method_train <- selfTrainingG

  }

  else if (method %in% c("COBCG", "COBC")) {
    method_train <- coBCG
  }

  else if (method %in% c("DEMOCRATICG", "DEMOCRATIC")) {
    method_train <- democraticG
  }

  else if (method %in% c("SETREDG", "SETRED")) {
    method_train <- setredG
  }

  else if (method %in% c("SNNRCEG", "SNNRCE")) {
    method_train <- snnrceG
  }

  else if (method %in% c("TRITRAININGG", "TRITRAINING")) {
    method_train <- triTrainingG
  }

  method_train

}

#' @title FUNCTION TO TRAIN GENERIC MODEL
#' @description FUNCTION TO TRAIN GENERIC MODEL
#' @param y (optional) factor (classes)
#' @param ... list parms trControl (method...)
#' @return  model trained
#' @export
train_generic <- function(y, ...) {

  method <- trControl <- NULL

  #GET PARAMETERS IN LIST
  parms = list(...)
  #print(length(parms))
  #print(names(parms))

  #CHECK IF EXISTS PARAMETERS
  if (length(parms) > 0) {

    #ASSIGN BY NAME PARAM
    for (name in names(parms)) {
      assign(name, parms[[name]])
    }

    model_train <- NULL
    #GET FUNCTION GENERIC
    method_train <- get_function_generic(method)

    #CALL WITH PARAMETERS
    model_train <- do.call(method_train, c(list(y = y), trControl))
    model_train
  }

  #NOT EXISTS PARMS
  else {
    stop("Not exists parms")
  }
}





new_model_sslr <- function(train_function, model, args) {

  result <- list(
    fit_function = train_function,
    model = model,
    args = args
  )
  class(result) <- c("model_sslr")

  result
}

new_model_sslr_fitted <- function(model, model_name, args,colnames_x,elapsed, formula = NULL) {

  result <- list(
    model = model,
    model_name = model_name,
    args = args,
    col_names_fit = colnames_x,
    elapsed = elapsed,
    formula = formula
  )
  class(result) <- c("model_sslr_fitted")

  result
}


#' @title Print model SSLR
#' @param object model_sslr object to print
#' @importFrom purrr map
print.model_sslr <- function(object) {

  cat("Model is:", object$model)

  cat("\n\nArgs:")
  eng_args <- map(object$args, convert_arg)
  cat(print_arg_list(eng_args), "\n", sep = "")
}


convert_arg <- function(x) {
  x
}

deparserizer <- function(x, limit = options()$width - 10) {
  x <- deparse(x, width.cutoff = limit)
  x <- gsub("^    ", "", x)
  x <- paste0(x, collapse = "")
  if (nchar(x) > limit)
    x <- paste0(substring(x, first = 1, last = limit - 7), "<snip>")
  x
}

print_arg_list <- function(x, ...) {
  atomic <- vapply(x, is.atomic, logical(1))
  x2 <- x
  x2[!atomic] <- lapply(x2[!atomic], deparserizer, ...)
  res <- paste0("  ", names(x2), " = ", x2, collaspe = "\n")
  cat(res, sep = "")
}



#' @title Fit with x and y
#' @description Funtion to fit with x and y
#' @param object is the model
#' @param x is a data frame or matrix with train dataset without objective feature.
#' X have labeled and unlabeled data
#' @param y is objective feature with labeled values and NA values in unlabeled data
#' @param ... unused in this case
#' @export fit_xy.model_sslr
#' @export
fit_xy.model_sslr <- function(object, x = NULL, y = NULL, ...) {

  eval_env <- rlang::env()
  eval_env$x <- x
  eval_env$y <- y

  fit_interface <- check_xy_interface(eval_env$x, eval_env$y)

  elapsed <- system.time(model <- object$fit_function(eval_env$x, eval_env$y))

  new_model_sslr_fitted(model,class(model),object$args,colnames(eval_env$x),elapsed,NULL)

}



#' @title Fit with formula and data
#' @description Funtion to fit through the formula
#' @param object is the model
#' @param formula is the formula
#' @param data is the total data train
#' @param ... unused in this case
#' @importFrom rlang quos
#' @export fit.model_sslr
#' @export
fit.model_sslr <- function(object, formula = NULL, data = NULL, ...) {

  dots <- quos(...)


  if (all(c("x", "y") %in% names(dots)))
    rlang::abort("`fit.model_sslr()` is for the formula methods. Use `fit_xy()` instead.")

  fit_form_interface <- check_form_interface(formula, data)

  x_and_y <- get_x_y(formula, data)

  eval_env <- rlang::env()
  eval_env$x <- x_and_y$x
  eval_env$y <- x_and_y$y

  fit_interface <- check_xy_interface(eval_env$x, eval_env$y)

  elapsed <- system.time(model <- object$fit_function(eval_env$x, eval_env$y))

  new_model_sslr_fitted(model,class(model),object$args,colnames(eval_env$x),elapsed,formula)

}


#' fit_x_u
#' @title fit_x_u object
#' @param object object
#' @param ... other parameters to be passed
#' @export
fit_x_u <- function(object, ...){
  UseMethod("fit_x_u")
}


#' @title Fit with x , y (labeled data) and unlabeled data (x_U)
#' @description Funtion to fit with x and y and x_U.
#' Function calcule y with NA values and append in y param
#' @param object is the model
#' @param x is a data frame or matrix with train dataset without objective feature.
#' X only have labeled data
#' @param y is objective feature with labeled values
#' @param x_U train unlabeled data without objective feature
#' @param ... This parameter is included for compatibility reasons.
#' @export fit_x_u.model_sslr
#' @export
fit_x_u.model_sslr <- function(object, x = NULL, y = NULL, x_U = NULL, ...) {

  eval_env <- rlang::env()
  eval_env$x <- rbind(x, x_U)

  y_unlabeled <- rep(NA, nrow(x_U))

  if (is.factor(y))
    eval_env$y <- factor(append(as.character(y), y_unlabeled))

  else
    eval_env$y <- as.numeric(append(as.numeric(y), y_unlabeled))


  fit_interface <- check_xy_interface(eval_env$x, eval_env$y)

  elapsed <- system.time(model <- object$fit_function(eval_env$x, eval_env$y))

  new_model_sslr_fitted(model,class(model),object$args,colnames(eval_env$x),elapsed,NULL)

}




inher <- function(x, cls) {
  if (!is.null(x) && !inherits(x, cls)) {
    call <- match.call()
    obj <- deparse(call[["x"]])
    if (length(cls) > 1)
      stop(
        "`", obj, "` should be one of the following classes: ",
        paste0("'", cls, "'", collapse = ", "), call. = FALSE
      )
    else
      stop(
        "`", obj, "` should be a ", cls, " object", call. = FALSE
      )
  }
  invisible(x)
}


#' @title Ceck interface x y
#' @description Check interface
#' @param x data without class labels
#' @param y values class
check_xy_interface <- function(x, y) {

  # `y` can be a vector (which is not a class), or a factor (which is not a vector)
  if (!is.null(y) && !is.vector(y))
    inher(y, c("data.frame", "matrix", "factor", "numeric"))

  # Determine the `fit()` interface
  matrix_interface <- !is.null(x) & !is.null(y) && is.matrix(x)
  df_interface <- !is.null(x) & !is.null(y) && is.data.frame(x)

  if (nrow(x) != length(y)) {
    stop("Length x Instances is not equal length y")
  }

  #Check labeled and unlabeled data
  check_labeled_unlabeled(y)

  if (matrix_interface)
    return("data.frame")
  if (df_interface)
    return("data.frame")
  stop("Error when checking the interface")
}


check_form_interface <- function(formula, data) {
  inher(formula, "formula")
  inher(data, c("data.frame", "tbl_spark", "matrix"))

  # Determine the `fit()` interface
  form_interface <- !is.null(formula) & !is.null(data)

  if (form_interface)
    return("formula")
  stop("Error when checking the interface")
}


check_labeled_unlabeled <- function(y){
  # Obtain the indexes of labeled and unlabeled instances
  labeled <- which(!is.na(y))
  unlabeled <- which(is.na(y))
  ## Check the labeled and unlabeled sets
  if (length(labeled) == 0) {
    # labeled is empty
    stop("The labeled set is empty. All the values in y parameter are NA.")
  }
  if (length(unlabeled) == 0) {
    # unlabeled is empty
    stop("The unlabeled set is empty. None value in y parameter is NA.")
  }
}


get_doParallel_loaded <- function() {
  "doParallel" %in% (.packages())
}


get_x_y_And_unlabeled <- function(X, y) {

  index_labeled <- which(!is.na(y))

  list(x = X[index_labeled,], y = y[index_labeled], X_u = X[-index_labeled,])

}

#' @importFrom dplyr tibble
pred_factor_tibble <- function(x){
    tibble(.pred_class = x)
}

#' @importFrom dplyr tibble
pred_numeric_tibble <- function(x){
  tibble(.pred = x)
}

#' @title Load RSSL
#' @description function to load RSSL package
load_RSSL <- function(){

  cond <- "RSSL" %in% (.packages())

  if(!cond){
    if (requireNamespace("RSSL", quietly=TRUE)){

    }
    else{
      rlang::abort("You should install RSSL package with install.packages(RSSL)")
    }
  }

}


#' @title Load conclust
#' @description function to load conclust package
load_conclust <- function(){

  cond <- "conclust" %in% (.packages())

  if(!cond){
    if (requireNamespace("conclust", quietly=TRUE)){

    }
    else{
      rlang::abort("You should install conclust package with install.packages(conclust)")
    }
  }

}


#' @title Load parsnip
#' @description function to load parsnip package
load_parsnip <- function(){

  cond <- "parsnip" %in% (.packages())

  if(!cond){
    if (requireNamespace("parsnip", quietly=TRUE)){

    }
    else{
      rlang::abort("You should install parsnip package with install.packages(parsnip)
                   or install.packages(tidymodels)")
    }
  }

}


#' @title Load parsnip
#' @description function to load parsnip package
load_RANN <- function(){

  cond <- "RANN" %in% (.packages())

  if(!cond){
    if (requireNamespace("RANN", quietly=TRUE)){

    }
    else{
      rlang::abort("You should install RANN package with install.packages(RANN)")
    }
  }

}


#' Cluster labels
#' @title Get labels of clusters
#' @param object object
#' @param ... other parameters to be passed
#' @export
cluster_labels <- function(object, ...){
  UseMethod("cluster_labels")
}




#' @title Cluster labels
#' @description Get labels of clusters
#' raw returns factor or numeric values
#' @param object  model_sslr_fitted model built
#' @param type of predict in principal model: class, raw
#' @param ... other parameters to be passed
#' @export cluster_labels.model_sslr_fitted
#' @export
cluster_labels.model_sslr_fitted <- function(object,type = "class",...){

  if(object$model$mode != "clustering")
    stop("This model is not a clustering model")

  result <- as.factor(object$model$cluster)

  if(is.factor(result) & type == "class")
    result <- pred_factor_tibble(result)

  result


}


#' Centers clustering
#' @title Get centers model of clustering
#' @param object object
#' @param ... other parameters to be passed
#' @export
get_centers <- function(object, ...){
  UseMethod("get_centers")
}


#' @title Cluster labels
#' @description Get labels of clusters
#' raw returns factor or numeric values
#' @param object  model_sslr_fitted model built
#' @param ... other parameters to be passed
#' @export get_centers.model_sslr_fitted
#' @export
get_centers.model_sslr_fitted <- function(object, ...){

  if(object$model$mode != "clustering")
    stop("This model is not a clustering model")

  object$model$centers


}
