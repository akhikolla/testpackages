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
#  c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE,
#                                   seed = 25)
#  
#  # Train model
#  model <- train("Maxent", data = train)
#  
#  # Train cross validation model
#  folds <- randomFolds(data, k = 4, only_presence = TRUE, seed = 25)
#  cv_model <- train("Maxent", data = data, folds = folds)

## ----maxent results-----------------------------------------------------------
#  model@model@results

## ----maxent var importance----------------------------------------------------
#  vi <- maxentVarImp(model)
#  vi

## ----maxent percent contribution plot-----------------------------------------
#  plotVarImp(vi[, 1:2])

## ----maxent permutation importance plot---------------------------------------
#  plotVarImp(vi[, c(1,3)])

## ----maxnet model-------------------------------------------------------------
#  maxnet_model <- train("Maxnet", data = train)

## -----------------------------------------------------------------------------
#  vi_maxnet <- varImp(maxnet_model, permut = 5)
#  vi_maxnet

## ----plot var imp-------------------------------------------------------------
#  plotVarImp(vi_maxnet)

## ----permut exercise----------------------------------------------------------
#  # Compute the permutation importance
#  vi_maxent <- varImp(model, permut = 10)
#  # Print it
#  vi_maxent
#  # Compare with Maxent output
#  maxentVarImp(model)

## ----jk-----------------------------------------------------------------------
#  jk <- doJk(maxnet_model, metric = "auc", test = test)
#  jk

## ----plot jk train------------------------------------------------------------
#  plotJk(jk, type = "train", ref = auc(maxnet_model))

## ----plot jk test, fig.align="center"-----------------------------------------
#  plotJk(jk, type = "test", ref = auc(maxnet_model, test = test))

## ----plot bio1, fig.align="center"--------------------------------------------
#  plotResponse(maxnet_model, var = "bio1", type = "cloglog", only_presence = TRUE,
#               marginal = FALSE, rug = TRUE)

## ----plot biome, fig.align="center"-------------------------------------------
#  plotResponse(maxnet_model, var = "biome", type = "logistic",
#               only_presence = TRUE, marginal = TRUE, fun = mean, color = "blue")

## ----plot cv response---------------------------------------------------------
#  plotResponse(cv_model, var = "bio1", type = "cloglog", only_presence = TRUE,
#               marginal = TRUE, fun = mean, rug = TRUE)

## ----report-------------------------------------------------------------------
#  modelReport(maxnet_model, type = "cloglog", folder = "virtual-sp", test = test,
#              response_curves = TRUE, only_presence = TRUE, jk = TRUE,
#              env = predictors)

## ----backgrounds--------------------------------------------------------------
#  set.seed(25)
#  bg <- dismo::randomPoints(predictors, 10000)
#  bg <- prepareSWD(species = "Bgs", a = bg, env = predictors,
#                   categorical = "biome")

## ----heat map-----------------------------------------------------------------
#  plotCor(bg, method = "spearman", cor_th = 0.7)

## ----cor var------------------------------------------------------------------
#  corVar(bg, method = "spearman", cor_th = 0.7)

## ----varSel-------------------------------------------------------------------
#  selected_variables_model <- varSel(maxnet_model, metric = "auc", test = test,
#                                     bg4cor = bg, method = "spearman",
#                                     cor_th = 0.7, permut = 1)

## ----output varSel------------------------------------------------------------
#  selected_variables_model

## ----exercise 1---------------------------------------------------------------
#  selected_variables_model <- varSel(model, metric = "aicc", bg4cor = bg,
#                                     method = "spearman", cor_th = 0.7,
#                                     env = predictors, use_pc = TRUE)

## ----permutation--------------------------------------------------------------
#  varImp(model, permut = 1)

## ----reduce var 1-------------------------------------------------------------
#  cat("Testing AUC before: ", auc(maxnet_model, test = test))
#  reduced_variables_model <- reduceVar(maxnet_model, th = 6, metric = "auc",
#                                       test = test, permut = 1)
#  cat("Testing AUC after: ", auc(reduced_variables_model, test = test))

## ----reduce var 2-------------------------------------------------------------
#  cat("Testing AUC before: ", auc(maxnet_model, test = test))
#  reduced_variables_model <- reduceVar(maxnet_model, th = 15, metric = "auc",
#                                       test = test, permut = 1, use_jk = TRUE)
#  cat("Testing AUC after: ", auc(reduced_variables_model, test = test))

