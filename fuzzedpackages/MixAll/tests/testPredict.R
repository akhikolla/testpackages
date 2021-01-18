library(MixAll)

testPredict<-function(nbTrain , nbTest)
{
  ## test categorical predictions
  train1 <- matrix( c( sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.05,0.05,0.9))
                     , sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.9,0.05,0.05))
                     , sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.05,0.05,0.9))
                     , sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.9,0.05,0.05))
                     )
                  , ncol =2
                  )
  model <- clusterCategorical(train1,2,models = "categorical_p_pjk")
  pred  <- clusterPredict(train1,model)
  # more than 5 classification errors is abnormal
  if ( sum(pred@zi - model@zi) != 0)
  { print("Predict Categorical failed");return(FALSE)}
  
  ##------------------------------------------------------------------------------
  ## test Poisson predictions
  train2 <- matrix( c( rpois(nbTrain,lambda = 1), rpois(nbTrain,lambda = 10)
                     , rpois(nbTrain,lambda = 1), rpois(nbTrain,lambda = 10)
                     )
                  , ncol =2
                  )
  model <- clusterPoisson(train2,2,models = "poisson_p_lk")
  pred  <- clusterPredict(train2,model)
  # Predictions should be the same
  if ( sum(pred@zi -model@zi) != 0)
  { print("Predict Poisson failed");return(FALSE)}
  
  ##------------------------------------------------------------------------------
  ## test Gaussian predictions
  train3 <- matrix( c( rnorm(nbTrain, mean = 1, sd=1), rnorm(nbTrain,mean = 10, sd=1)
                     , rnorm(nbTrain, mean = 1, sd=1), rnorm(nbTrain,mean = 10, sd=1)
                     )
                  , ncol =2
                  )
  model <- clusterDiagGaussian(train3, 2, models = "gaussian_p_s")
  pred  <- clusterPredict(train3, model)
  # Predictions should be the same
  if ( sum(pred@zi -model@zi) != 0)
  { print("Predict Gaussian failed");return(FALSE)}
  
  ##------------------------------------------------------------------------------
  ## test gamma predictions
  train4 <- matrix( c( rgamma(nbTrain, shape = 1, scale=1), rgamma(nbTrain,shape = 10, scale=1)
                     , rgamma(nbTrain, shape = 1, scale=1), rgamma(nbTrain,shape = 10, scale=1)
                     )
                  , ncol =2
                  )
  model <- clusterGamma(train4, 2, models = "gamma_p_ak_b")
  pred  <- clusterPredict(train4,model)
  # more than 5 classification errors is abnormal
  if ( sum(pred@zi -model@zi) != 0)
  { print("Predict gamma failed");return(FALSE)}
  
  ##------------------------------------------------------------------------------
  ## test mixed data predictions
  train <- list(train1, train2, train3, train4)
  models <- c("categorical_p_pjk", "poisson_p_lk",  "gaussian_p_s","gamma_p_ak_b")
  
  model <- clusterMixedData(train, models, 2)
  pred  <- clusterPredict(train, model)
  # more than 5 classification errors is abnormal
  if ( sum(pred@zi -model@zi) != 0)
  { print("Predict mixed failed");return(FALSE)}
  
  ##------------------------------------------------------------------------------
  return(TRUE)  
}

testPredict(1000, 20)
