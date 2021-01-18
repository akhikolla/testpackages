#' @title Generic Interface coBCReg model
#' @description coBCReg is based on an
#' ensemble of N diverse regressors. At each iteration and for each regressor, the
#' companion committee labels the unlabeled examples then the regressor select
#' the most informative newly-labeled examples for itself, where the selection confidence
#' is based on estimating the validation error. The final prediction is the
#' average of the estimates of the N regressors.
#' @details For regression tasks, labeling data is very expensive computationally. Its so slow.
#' @param y A vector with the labels of training instances. In this vector the
#' unlabeled instances are specified with the value \code{NA}.
#' @param gen.learner A function for training \code{N} supervised base classifiers.
#' This function needs two parameters, indexes and cls, where indexes indicates
#' the instances to use and cls specifies the classes of those instances.
#' @param gen.pred A function for predicting the probabilities per classes.
#' This function must be two parameters, model and indexes, where the model
#' is a classifier trained with \code{gen.learner} function and
#' indexes indicates the instances to predict.
#' @param N The number of classifiers used as committee members. All these classifiers
#' are trained using the \code{gen.learner} function. Default is 3.
#' @param perc.full A number between 0 and 1. If the percentage
#' of new labeled examples reaches this value the self-labeling process is stopped.
#' Default is 0.7.
#' @param u Number of unlabeled instances in the pool. Default is 100.
#' @param max.iter Maximum number of iterations to execute in the self-labeling process.
#' Default is 50.
#' @param gr growing rate
#' @references
#' Mohamed Farouk Abdel-Hady, Mohamed Farouk Abdel-Hady and Günther Palm.\cr
#' \emph{Semi-supervised Learning for Regression with Cotraining by Committee}\cr
#' Institute of Neural Information Processing
#' University of Ulm
#' D-89069 Ulm, Germany
coBCRegG <- function(
  y,
  gen.learner,
  gen.pred,
  N = 3,
  perc.full = 0.7,
  u = 100,
  max.iter = 50,
  gr = 1
) {

  # Check N
  if (N <= 0) {
    stop("Parameter N must be a positive number.")
  }
  # Check max.iter
  if (max.iter < 1) {
    stop("Parameter max.iter is less than 1. Expected a value greater than and equal to 1.")
  }
  # Check perc.full
  if (perc.full < 0 || perc.full > 1) {
    stop("Parameter perc.full is not in the range 0 to 1.")
  }

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

  #Init variables
  His<- vector(mode = "list", length = N)

  Lis<- vector(mode = "list", length = N)
  Vis <- vector(mode = "list", length = N)


  # Lists to store the training instances indexes and it's predicted values
  Lind <- vector(mode = "list", length = N)
  Lcls <- vector(mode = "list", length = N)

  #Create boostrap samples

  for(i in 1:N)
  {
    #BOOSTRAP SAMPLE
    Li <- sample(x = labeled, size = length(labeled), replace = TRUE)
    Vi <- setdiff(labeled, Li)

    Lis[[i]] <- Li
    Vis[[i]] <- Vi

    #Update List indexes and cls, cls is y and update each iteration
    Lind[[i]] <- c()
    Lcls[[i]] <- y

    #Train model
    Hi <- gen.learner(Li, y[Li])
    His[[i]] <- Hi

  }

  #Init iter
  iter <- 1
  min.amount <- round(length(unlabeled) * (1 - perc.full))

  #Variable to detect some model
  change_model <- FALSE

  while ((length(unlabeled) > min.amount) && (iter <= max.iter)) {

    #If unlabeled is empty
    if(length(unlabeled) == 0) break

    #First change model to FALSE
    change_model <- FALSE


    for(j in 1:N){

      # Select randomly a pool of unlabeled instances
      pool <- sample(x = unlabeled, size = min(u, length(unlabeled)), replace = FALSE)

      #Select relevant examples
      sel <- selectRelevantExamples(His,j,pool, Vis[[j]], gr, Lis[[j]], gen.learner,gen.pred,Lcls[[j]])

      #Delete select examples in unlabeled set
      unlabeled <- setdiff(unlabeled,sel)

      #Update indexes and predicted values (pi)
      Lind[[j]] <- sel

      #If selected values, update y in each model
      if(length(sel) > 0){
        values <- Lcls[[j]]
        values[sel] <- as.numeric((as.matrix(gen.pred(His[[j]], sel))[,1]))

        Lcls[[j]] <- values

        #Change model (selected examples > 0)
        change_model <- TRUE
      }

      else{
        Lind[[j]] <- c()
      }


    }


    for(j in 1:N){

      if(length(Lind[[j]]) > 0){

        #Add new labeled data
        Lis[[j]] <- unique(c(Lis[[j]],Lind[[j]]))

        #Train model
        current.y <- Lcls[[j]]

        Hi <- gen.learner(Lis[[j]], current.y[Lis[[j]]])
        His[[j]] <- Hi
      }

    }

    #If not detect any change, exit
    if(!change_model) break

    iter <- iter + 1

  }

  # Save result
  result <- list(
    model = His,
    gen.pred = gen.pred,
    model.index = Lis
  )
  class(result) <- "coBCRegG"

  return(result)

}



selectRelevantExamples <- function(models,pos,pool,Vi,gr,Li,gen.learner,gen.pred,y){

  #Calculate validation error using Vi

  model <- models[[pos]]

  predictions <- as.numeric((as.matrix(gen.pred(model,Vi))[,1]))

  e_j <- sqrt(mean((predictions - y[Vi])^2))

  #Create deltas vector (zeros)
  deltas <- rep(0, length(pool))

  for(i in 1:length(pool)){

    xu <- pool[i]

    prediction <- predict_xu(models,xu,gen.pred,pos)

    temp.L <- c(Li,xu)
    temp.Lcls <- y

    temp.Lcls[xu] <- prediction

    temp.h <- gen.learner(temp.L,temp.Lcls[temp.L])

    #Calculate e (validation error)

    predictions <- as.numeric((as.matrix(gen.pred(temp.h,Vi))[,1]))

    e <- sqrt(mean((predictions - temp.Lcls[Vi])^2))


    deltas[i] <- (e_j - e)/e_j

  }

  pi_temp <- c()

  for(i in 1:gr){

    #Get max delta
    idx_max_delta <-  which.max(deltas)

    if(length(idx_max_delta) == 0 ) break


    #If > 0
    if(deltas[idx_max_delta] > 0){

      pi_temp <- c(pi_temp,pool[idx_max_delta])

      #Update delta to -1 (under 0), not select again
      deltas[idx_max_delta] <- 0 #Zero, not select again
    }


  }


  return(pi_temp)

}


predict_xu <- function(models,xu,gen.pred, pos = 1){

  #Get num models
  n <- length(models)

  #If n = 1, predicted values if in first model
  if(n == 1){
    return(as.numeric((as.matrix(gen.pred(models[[1]],xu))[,1])))
  }

  #Predicted acumulator
  Hj <- 0

  for(j in 1:n){

    if(j == pos) next

    Hj <- Hj + as.numeric((as.matrix(gen.pred(models[[j]],xu))[,1]))

  }

  #Final value
  return(Hj/(n-1))

}



#' @title General Interface coBCReg model
#' @description coBCReg is based on an
#' ensemble of N diverse regressors. At each iteration and for each regressor, the
#' companion committee labels the unlabeled examples then the regressor select
#' the most informative newly-labeled examples for itself, where the selection confidence
#' is based on estimating the validation error. The final prediction is the
#' average of the estimates of the N regressors.
#' @details For regression tasks, labeling data is very expensive computationally. Its so slow.
#' @param learner model from parsnip package for training a supervised base classifier
#' using a set of instances. This model need to have probability predictions
#' @param N The number of classifiers used as committee members. All these classifiers
#' are trained using the \code{gen.learner} function. Default is 3.
#' @param perc.full A number between 0 and 1. If the percentage
#' of new labeled examples reaches this value the self-labeling process is stopped.
#' Default is 0.7.
#' @param u Number of unlabeled instances in the pool. Default is 100.
#' @param max.iter Maximum number of iterations to execute in the self-labeling process.
#' Default is 50.
#' @references
#' Mohamed Farouk Abdel-Hady, Mohamed Farouk Abdel-Hady and Günther Palm.\cr
#' \emph{Semi-supervised Learning for Regression with Cotraining by Committee}\cr
#' Institute of Neural Information Processing
#' University of Ulm
#' D-89069 Ulm, Germany
coBCReg <- function(
  learner,
  N = 3,
  perc.full = 0.7,
  u = 100,
  max.iter = 50
) {

  train_function <- function(x, y) {
    ### Check parameters ###
    x <- as.data.frame(x)
    y <- as.numeric(y)

    # Instance matrix case
    gen.learner1 <- function(training.ints, cls) {
      m <- learner %>% parsnip::fit_xy(x = x[training.ints,], y = cls)
      return(m)
    }

    gen.pred1 <- function(m, testing.ints) {
      prob <- predict(m, x[testing.ints,], type = "numeric")
      return(prob)
    }

    result <- coBCRegG(y, gen.learner1, gen.pred1, N, perc.full, u, max.iter)

    ### Result ###
    result$pred.params = c("numeric","raw")
    result$mode = "regression"
    class(result) <- "coBCReg"

    return(result)
  }

  args <- list(
    learner = learner,
    N = N,
    perc.full = perc.full,
    u = u,
    max.iter = max.iter
  )

  new_model_sslr(train_function, "coBCReg", args)
}



#' @method predict coBCReg
#' @importFrom stats predict
predict.coBCReg <- function(object, x,...) {

  His <- object$model
  N <- length(His)

  nrow_data <- nrow(x)
  results <- matrix(data = rep(0,nrow_data*N), nrow =  nrow_data)


  for(i in 1:N){

    model <- His[[i]]

    results[,i] <- as.numeric((as.matrix(predict(model, x,"numeric"))[,1]))

  }

  rowMeans(results)

}
