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
#  p_coords <- virtualSp$presence
#  a_coords <- virtualSp$absence
#  data <- prepareSWD(species = "Virtual species", p = p_coords, a = a_coords,
#                     env = predictors[[1:8]])
#  
#  # Split data in training and testing datasets
#  c(train, test) %<-% trainValTest(data, test = 0.2, seed = 25)
#  cat("# Training  : ", nrow(train@data))
#  cat("\n# Testing   : ", nrow(test@data))
#  
#  # Create folds
#  folds <- randomFolds(train, k = 4, seed = 25)

## -----------------------------------------------------------------------------
#  set.seed(25)
#  model <- train("ANN", data = train, size = 10, folds = folds)
#  model

## ----auc----------------------------------------------------------------------
#  auc(model)
#  auc(model, test = TRUE)

## ----get tunable args---------------------------------------------------------
#  getTunableArgs(model)

## ----optimize model-----------------------------------------------------------
#  h <- list(size = 10:50, decay = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
#            maxit = c(50, 100, 300, 500))
#  
#  om <- optimizeModel(model, hypers = h, metric = "auc", seed = 25)

## ----best model---------------------------------------------------------------
#  best_model <- om@models[[1]]
#  om@results[1, ]

## ----evaluate final model, fig.align="center"---------------------------------
#  set.seed(25)
#  final_model <- train("ANN", data = train, size = om@results[1, 1],
#                       decay = om@results[1, 2], maxit = om@results[1, 4])
#  plotROC(final_model, test = test)

