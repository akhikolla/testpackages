## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE, eval = FALSE,
                      fig.align = "center")

## ----load data----------------------------------------------------------------
#  library(SDMtune)
#  library(zeallot)
#  
#  # Prepare data
#  files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
#                      pattern = "grd", full.names = TRUE)
#  predictors <- raster::stack(files)
#  data <- prepareSWD(species = "Virtual species", p = virtualSp$presence,
#                     a = virtualSp$background, env = predictors,
#                     categorical = "biome")
#  
#  # Split data in training, validation and testing datasets
#  c(train, val, test) %<-% trainValTest(data, val = 0.2, test = 0.2,
#                                        only_presence = TRUE, seed = 61516)
#  cat("# Training  : ", nrow(train@data))
#  cat("\n# Validation: ", nrow(val@data))
#  cat("\n# Testing   : ", nrow(test@data))
#  
#  # Train Maxnet model with default settings
#  model <- train("Maxnet", data = train)

## ----tune bg------------------------------------------------------------------
#  # Define the values for bg
#  h <- list(reg = seq(0.2, 1, 0.1))
#  # Call the gridSearch function
#  exp_1 <- gridSearch(model, hypers = h, metric = "auc", test = val)

## ----print exp_1--------------------------------------------------------------
#  exp_1

## ----plot exp_1---------------------------------------------------------------
#  plot(exp_1, title = "Experiment 1")

## ----plot exp_1 iter----------------------------------------------------------
#  plot(exp_1, title = "Experiment 1", interactive = TRUE)

## ----slot results-------------------------------------------------------------
#  exp_1@results

## ----order results------------------------------------------------------------
#  exp_1@results[order(-exp_1@results$test_AUC), ]

## ----exercise 1, eval=FALSE, class.source='exercise'--------------------------
#  # Define the values for reg
#  h <- list(reg = 1:4)
#  # Call the gridSearch function
#  exp_2 <- gridSearch(model, hypers = h, metric = "tss", test = val)

## ----exercise 2, eval=FALSE, class.source='exercise'--------------------------
#  # Define the values for fc
#  h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"))
#  # Call the gridSearch function
#  exp_3 <- gridSearch(model, hypers = h, metric = "auc", test = val)

## ----exercise 3, eval=FALSE, class.source='exercise'--------------------------
#  maxent_model <- train("Maxent", data = data)
#  # Define the values for fc
#  h <- list("iter" = seq(300, 1100, 200))
#  # Call the gridSearch function
#  exp_4 <- gridSearch(maxent_model, hypers = h, metric = "auc", test = val)

## ----get hypers---------------------------------------------------------------
#  getTunableArgs(model)

## ----exp 5, eval=FALSE--------------------------------------------------------
#  h <- list(reg = seq(0.2, 2, 0.2),
#            fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"))
#  exp_5 <- gridSearch(model, hypers = h, metric = "auc", test = val)

## ----random search------------------------------------------------------------
#  h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lp", "lqp", "lqph"))
#  exp_6 <- randomSearch(model, hypers = h, metric = "auc", test = val, pop = 10,
#                        seed = 65466)

## ----sdmtune results----------------------------------------------------------
#  exp_6@results

## ----optimize model, eval=FALSE-----------------------------------------------
#  exp_7 <- optimizeModel(model, hypers = h, metric = "auc", test = val, pop = 15,
#                         gen = 2, seed = 798)

## ----merge--------------------------------------------------------------------
#  index <- which.max(exp_6@results$test_AUC)  # Index of the best model in the experiment
#  new_train <- exp_6@models[[index]]@data  # New train dataset containing only the selected variables
#  merged_data <- mergeSWD(new_train, val, only_presence = TRUE) # Merge only presence data

## ----final model--------------------------------------------------------------
#  final_model <- train("Maxnet", data = merged_data, fc = exp_6@results[index, 1],
#                       reg = exp_6@results[index, 2])

## ----final evaluation---------------------------------------------------------
#  auc(final_model, test = test)

## ----cross validation---------------------------------------------------------
#  # Create the folds from the training dataset
#  folds <- randomFolds(train, k = 4, only_presence = TRUE, seed = 25)
#  # Train the model
#  cv_model <- train("Maxent", data = train, folds = folds)

## ----randomSearch cv----------------------------------------------------------
#  h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lp", "lqp", "lqph"))
#  exp_8 <- randomSearch(cv_model, hypers = h, metric = "auc", pop = 10,
#                        seed = 65466)

## ----final evaluation cv------------------------------------------------------
#  final_model <- train("Maxent", data = exp_8@models[[1]]@data,
#                       fc = exp_8@results[1, 1], reg = exp_8@results[1, 2])
#  auc(final_model, test = test)

