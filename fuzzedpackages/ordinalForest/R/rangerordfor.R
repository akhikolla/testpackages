##' @author Marvin N. Wright
##' @useDynLib ordinalForest
##' @importFrom Rcpp evalCpp
##' @import stats 
##' @import utils
rangerordfor <- function(formula = NULL, data = NULL, num.trees = 500, mtry = NULL,
                   importance = "none", write.forest = TRUE, probability = FALSE,
                   min.node.size = NULL, replace = TRUE, 
                   sample.fraction = ifelse(replace, 1, 0.632), 
                   case.weights = NULL, 
                   splitrule = NULL, alpha = 0.5, minprop = 0.1,
                   split.select.weights = NULL, always.split.variables = NULL,
                   respect.unordered.factors = "ignore",
                   scale.permutation.importance = FALSE,
                   keep.inbag = FALSE, holdout = FALSE,
                   num.threads = NULL, save.memory = FALSE,
                   verbose = TRUE, seed = NULL, 
                   dependent.variable.name = NULL, status.variable.name = NULL, 
                   classification = NULL, borders = NULL, userps=FALSE) { 
  
  ## GenABEL GWA data
  if ("gwaa.data" %in% class(data)) {
    snp.names <- data@gtdata@snpnames
    sparse.data <- data@gtdata@gtps@.Data
    data <- data@phdata
    if ("id" %in% names(data)) {
      data$"id" <- NULL
    }
    gwa.mode <- TRUE
    save.memory <- FALSE
  } else {
    sparse.data <- as.matrix(0)
    gwa.mode <- FALSE
  }
  
  ## Check missing values
  if (any(is.na(data))) {
    offending_columns <- colnames(data)[colSums(is.na(data)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
  ## Formula interface. Use whole data frame is no formula provided and depvarname given
  if (is.null(formula)) {
    if (is.null(dependent.variable.name)) {
      stop("Error: Please give formula or dependent variable name.")
    }
    if (is.null(status.variable.name)) {
      status.variable.name <- "none"
      response <- data[, dependent.variable.name]
    } else {
      response <- data[, c(dependent.variable.name, status.variable.name)]
    }
    data.selected <- data
  } else {
    formula <- formula(formula)
    if (class(formula) != "formula") {
      stop("Error: Invalid formula.")
    }
    data.selected <- model.frame(formula, data, na.action = na.fail)
    response <- data.selected[[1]]
  }
  
  ## Treetype
  if (is.factor(response)) {
    if (probability) {
      treetype <- 9
    } else {
      treetype <- 1
    }
  } else if (is.numeric(response) & is.vector(response)) {
    if (!is.null(classification) && classification) {
      treetype <- 1
    } else if (probability) {
      treetype <- 9
    } else {
      treetype <- 3
    }
  } else if (class(response) == "Surv" | is.data.frame(response) | is.matrix(response)) {
    treetype <- 5
  } else {
    stop("Error: Unsupported type of dependent variable.")
  }
  
  ## Dependent and status variable name. For non-survival dummy status variable name.
  if (!is.null(formula)) {
    if (treetype == 5) {
      dependent.variable.name <- dimnames(response)[[2]][1]
      status.variable.name <- dimnames(response)[[2]][2]
    } else {
      dependent.variable.name <- names(data.selected)[1]
      status.variable.name <- "none"
    }
    independent.variable.names <- names(data.selected)[-1]
  } else {
    independent.variable.names <- colnames(data.selected)[colnames(data.selected) != dependent.variable.name &
                                                          colnames(data.selected) != status.variable.name]
  }
  
  ## Old version of if respect.unordered.factors
  if (respect.unordered.factors == TRUE) {
    respect.unordered.factors <- "order"
  } else if (respect.unordered.factors == FALSE) {
    respect.unordered.factors <- "ignore"
  }
  
  ## Recode characters as factors and recode factors if 'order' mode
  if (!is.matrix(data.selected)) {
    character.idx <- sapply(data.selected, is.character)
    
    if (respect.unordered.factors == "order") {
      ## Recode characters and unordered factors
      names.selected <- names(data.selected)
      ordered.idx <- sapply(data.selected, is.ordered)
      factor.idx <- sapply(data.selected, is.factor)
      independent.idx <- names.selected != dependent.variable.name & 
        names.selected != status.variable.name & 
        names.selected != paste0("Surv(", dependent.variable.name, ", ", status.variable.name, ")")
      recode.idx <- independent.idx & (character.idx | (factor.idx & !ordered.idx))

      ## Numeric response
      if (is.factor(response)) {
        num.response <- as.numeric(response)
      } else if (!is.null(dim(response))) {
        num.response <- response[, 1]
      } else {
        num.response <- response
      }

      ## Recode each column
      data.selected[recode.idx] <- lapply(data.selected[recode.idx], function(x) {
        ## Order factor levels
        means <- aggregate(num.response~x, FUN=mean)
        levels.ordered <- means$x[order(means$num.response)]
        
        ## Return reordered factor
        factor(x, levels = levels.ordered)
      })
      
      ## Save levels
      covariate.levels <- lapply(data.selected[independent.idx], levels)
    } else {
      ## Recode characters only
      data.selected[character.idx] <- lapply(data.selected[character.idx], factor)
    }
  }
  
  ## Input data and variable names, create final data matrix
  if (!is.null(formula) & treetype == 5) {
    data.final <- data.matrix(cbind(response[, 1], response[, 2],
                              data.selected[-1]))
    colnames(data.final) <- c(dependent.variable.name, status.variable.name,
                              independent.variable.names)
  } else if (is.matrix(data.selected)) {
    data.final <- data.selected
  } else {
    data.final <- data.matrix(data.selected)
  }
  variable.names <- colnames(data.final)
  
  ## If gwa mode, add snp variable names
  if (gwa.mode) {
    variable.names <- c(variable.names, snp.names)
    all.independent.variable.names <- c(independent.variable.names, snp.names)
  } else {
    all.independent.variable.names <- independent.variable.names
  }
  
  ## Number of trees
  if (!is.numeric(num.trees) | num.trees < 1) {
    stop("Error: Invalid value for num.trees.")
  }
  
  ## mtry
  if (is.null(mtry)) {
    mtry <- 0
  } else if (!is.numeric(mtry) | mtry < 0) {
    stop("Error: Invalid value for mtry")
  }
  
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  ## Keep inbag
  if (!is.logical(keep.inbag)) {
    stop("Error: Invalid value for keep.inbag")
  }
  
  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) | num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  
  ## Minumum node size
  if (is.null(min.node.size)) {
    min.node.size <- 0
  } else if (!is.numeric(min.node.size) | min.node.size < 0) {
    stop("Error: Invalid value for min.node.size")
  }
  
  ## Sample fraction
  if (!is.numeric(sample.fraction) | sample.fraction <= 0 | sample.fraction > 1) {
    stop("Error: Invalid value for sample.fraction. Please give a value in (0,1].")
  }
  
  ## Importance mode
  if (is.null(importance) | importance == "none") {
    importance.mode <- 0
  } else if (importance == "impurity") {
    importance.mode <- 1
    if (treetype == 5) {
      stop("Node impurity variable importance not supported for survival forests.")
    }
  } else if (importance == "permutation") {
    if (scale.permutation.importance) {
      importance.mode <- 2
    } else {
      importance.mode <- 3
    }
  } else {
    stop("Error: Unknown importance mode.")
  }
  
  ### ROMANNEU
  if(is.null(borders))
    borders <- c(0,0)
  
  ## Case weights: NULL for no weights
  if (is.null(case.weights)) {
    case.weights <- list(c(0,0)) # Kommentar Roman, zuvor: case.weights <- c(0,0)
    use.case.weights <- FALSE
    if (holdout) {
      stop("Error: Case weights required to use holdout mode.")
    }  # Kommentar Roman:
  } else if (is.numeric(case.weights)) {
    ## Sample from non-zero weights in holdout mode
    if (holdout) {
      sample.fraction <- sample.fraction * mean(case.weights > 0)
    }
    
    if (!replace && sum(case.weights > 0) < sample.fraction * nrow(data.final)) {
      stop("Error: Fewer non-zero case weights than observations to sample.")
    }
    case.weights <- list(case.weights)
    use.case.weights <- TRUE
  } else if (is.list(case.weights)) {
    if (length(case.weights) != num.trees) {
      stop("Error: Size of case weights list not equal to number of trees.")
    }
	## Sample from non-zero weights in holdout mode
    if (holdout) {
      warning("Warning: holdout set to FALSE, because there are tree-specific case weights.")
	  holdout <- FALSE
    }
    use.case.weights <- TRUE
  } else {
    stop("Error: Invalid case weights.")
  }
  
  if (use.case.weights & !replace) {
    if (any(sapply(case.weights, function(x) {mean(x == 0) < sample.fraction}))) {
      stop("Error: Too many 0's in case weights to draw without replacement.")
    }
  }

  # Kommentar Roman:
  # } else {
  #   use.case.weights <- TRUE
  #   
  #   ## Sample from non-zero weights in holdout mode
  #   if (holdout) {
  #     sample.fraction <- sample.fraction * mean(case.weights > 0)
  #   }
  #   
  #   if (!replace && sum(case.weights > 0) < sample.fraction * nrow(data.final)) {
  #     stop("Error: Fewer non-zero case weights than observations to sample.")
  #   }
  # }
  
  
  ## Split select weights: NULL for no weights
  if (is.null(split.select.weights)) {
    split.select.weights <- list(c(0,0))
    use.split.select.weights <- FALSE
  } else if (is.numeric(split.select.weights)) {
    if (length(split.select.weights) != length(all.independent.variable.names)) {
      stop("Error: Number of split select weights not equal to number of independent variables.")
    }
    split.select.weights <- list(split.select.weights)
    use.split.select.weights <- TRUE
  } else if (is.list(split.select.weights)) {
    if (length(split.select.weights) != num.trees) {
      stop("Error: Size of split select weights list not equal to number of trees.")
    }
    use.split.select.weights <- TRUE
  } else {
    stop("Error: Invalid split select weights.")
  }
  
  ## Always split variables: NULL for no variables
  if (is.null(always.split.variables)) {
    always.split.variables <- c("0", "0")
    use.always.split.variables <- FALSE
  } else {
    use.always.split.variables <- TRUE
  }
  
  if (use.split.select.weights & use.always.split.variables) {
    stop("Error: Please use only one option of split.select.weights and always.split.variables.")
  }
  
  ## Splitting rule
  if (is.null(splitrule)) {
    if (treetype == 5) {
      splitrule <- "logrank"
    } else if (treetype == 3) {
      splitrule <- "variance"
    }
    splitrule.num <- 1
  } else if (splitrule == "logrank") {
    if (treetype == 5) {
      splitrule.num <- 1
    } else {
      stop("Error: logrank splitrule applicable to survival data only.")
    }
  } else if (splitrule == "variance") {
    if (treetype == 3) {
      splitrule.num <- 2
    } else {
      stop("Error: variance splitrule applicable to regression data only.")
    }
  } else if (splitrule == "auc" | splitrule == "C") {
    if (treetype == 5) {
      splitrule.num <- 2
    } else {
      stop("Error: C index splitrule applicable to survival data only.")
    }
  } else if (splitrule == "auc_ignore_ties" | splitrule == "C_ignore_ties") {
    if (treetype == 5) {
      splitrule.num <- 3
    } else {
      stop("Error: C index splitrule applicable to survival data only.")
    }
  } else if (splitrule == "maxstat") {
    if (treetype == 5 | treetype == 3) {
      splitrule.num <- 4
    } else {
      stop("Error: maxstat splitrule applicable to regression or survival data only.")
    }
  } else {
    stop("Error: Unknown splitrule.")
  }
  
  ## Maxstat splitting
  if (alpha < 0 | alpha > 1) {
    stop("Error: Invalid value for alpha, please give a value between 0 and 1.")
  }
  if (minprop < 0 | minprop > 0.5) {
    stop("Error: Invalid value for minprop, please give a value between 0 and 0.5.")
  }

  ## Unordered factors  
  if (respect.unordered.factors == "partition") {
    names.selected <- names(data.selected)
    ordered.idx <- sapply(data.selected, is.ordered)
    factor.idx <- sapply(data.selected, is.factor)
    independent.idx <- names.selected != dependent.variable.name & names.selected != status.variable.name
    unordered.factor.variables <- names.selected[factor.idx & !ordered.idx & independent.idx]
    
    if (length(unordered.factor.variables) > 0) {
      use.unordered.factor.variables <- TRUE
      ## Check level count
      num.levels <- sapply(data.selected[, factor.idx & !ordered.idx & independent.idx, drop = FALSE], nlevels)
      max.level.count <- 8*.Machine$sizeof.pointer - 1
      if (max(num.levels) > max.level.count) {
        stop(paste("Too many levels in unordered categorical variable ", unordered.factor.variables[which.max(num.levels)], 
                   ". Only ", max.level.count, " levels allowed on this system. Consider using the 'order' option.", sep = ""))
      } 
    } else {
      unordered.factor.variables <- c("0", "0")
      use.unordered.factor.variables <- FALSE
    } 
  } else if (respect.unordered.factors == "ignore" | respect.unordered.factors == "order") {
    ## Ordering for "order" is handled above
    unordered.factor.variables <- c("0", "0")
    use.unordered.factor.variables <- FALSE
  } else {
    stop("Error: Invalid value for respect.unordered.factors, please use 'order', 'partition' or 'ignore'.")
  }

  ## Unordered maxstat splitting not possible
  if (use.unordered.factor.variables & !is.null(splitrule)) {
    if (splitrule == "maxstat") {
      stop("Error: Unordered factor splitting not implemented for 'maxstat' splitting rule.")
    } else if (splitrule %in% c("C", "auc", "C_ignore_ties", "auc_ignore_ties")) {
      stop("Error: Unordered factor splitting not implemented for 'C' splitting rule.")
    }
  }
  
  ## Warning for experimental 'order' splitting 
  if (respect.unordered.factors == "order") {
    if (treetype == 5) {
      warning("Warning: The 'order' mode for unordered factor handling for survival outcomes is experimental.")
    } else if (treetype == 1 | treetype == 9) {
      if (nlevels(response) > 2) {
        warning("Warning: The 'order' mode for unordered factor handling for multiclass classification is experimental.")
      }
    } else if (treetype == 3 & splitrule == "maxstat") {
      warning("Warning: The 'order' mode for unordered factor handling with the 'maxstat' splitrule is experimental.")
    }
  }

  ## Prediction mode always false. Use predict.ranger() method.
  prediction.mode <- FALSE
  predict.all <- FALSE
  prediction.type <- 1
  
  ## No loaded forest object
  loaded.forest <- list()
  
  ## Clean up
  rm("data.selected")
  
  ## Call Ranger
  result <- rangerCpp(treetype, dependent.variable.name, data.final, variable.names, mtry,
                      num.trees, verbose, seed, num.threads, write.forest, importance.mode,
                      min.node.size, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      status.variable.name, prediction.mode, loaded.forest, sparse.data,
                      replace, probability, unordered.factor.variables, use.unordered.factor.variables, 
                      save.memory, splitrule.num, case.weights, use.case.weights, predict.all, 
                      keep.inbag, sample.fraction, alpha, minprop, holdout, prediction.type, borders, userps)
  
  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Prepare results
  result$predictions <- drop(do.call(rbind, result$predictions))
  if (importance.mode != 0) {
    names(result$variable.importance) <- all.independent.variable.names
  }
  
  ## Set predictions
  if (treetype == 1 & is.factor(response)) {
    result$predictions <- integer.to.factor(result$predictions,
                                            levels(response))
    true.values <- integer.to.factor(unlist(data.final[, dependent.variable.name]),
                                     levels(response))
    result$confusion.matrix <- table(true.values, result$predictions, dnn = c("true", "predicted"))
  } else if (treetype == 5) {
    result$chf <- result$predictions
    result$predictions <- NULL
    result$survival <- exp(-result$chf)
  } else if (treetype == 9 & !is.matrix(data)) {
    ## Set colnames and sort by levels
    colnames(result$predictions) <- unique(response)
    result$predictions <- result$predictions[, levels(droplevels(response))]
  }
  
  ## Splitrule
  if (treetype == 3 | treetype == 5) {
    result$splitrule <- splitrule
  }
  
  ## Set treetype
  if (treetype == 1) {
    result$treetype <- "Classification"
  } else if (treetype == 3) {
    result$treetype <- "Regression"
  } else if (treetype == 5) {
    result$treetype <- "Survival"
  } else if (treetype == 9) {
    result$treetype <- "Probability estimation"
  }
  if (treetype == 3) {
    result$r.squared <- 1 - result$prediction.error / var(response)
  }
  result$call <- sys.call()
  result$importance.mode <- importance
  result$num.samples <- nrow(data.final)
  
  ## Write forest object
  if (write.forest) {
    result$forest$levels <- levels(response)
    result$forest$independent.variable.names <- independent.variable.names
    result$forest$treetype <- result$treetype
    class(result$forest) <- "ranger.forest"
    
    ## In 'ordered' mode, save covariate levels
    if (respect.unordered.factors == "order" & !is.matrix(data)) {
      result$forest$covariate.levels <- covariate.levels
    }
  }
  
  class(result) <- "ranger"
  return(result)
}




integer.to.factor <- function(x, labels) {
  factor(x, levels = seq_along(labels), labels = labels)
}
