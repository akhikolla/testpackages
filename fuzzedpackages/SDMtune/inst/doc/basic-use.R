## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE, eval = FALSE,
                      fig.align = "center")

## ----load pkgs----------------------------------------------------------------
#  library(ggplot2)    # To plot locations
#  library(maps)       # To access useful maps
#  library(rasterVis)  # To plot raster objects

## ----get predictors-----------------------------------------------------------
#  files <- list.files(path = file.path(system.file(package = "dismo"), "ex"),
#                      pattern = "grd", full.names = TRUE)

## ----create raster stack------------------------------------------------------
#  predictors <- raster::stack(files)

## ----names predicrtors--------------------------------------------------------
#  names(predictors)

## ----plot bio1----------------------------------------------------------------
#  gplot(predictors$bio1) +
#      geom_tile(aes(fill = value)) +
#      coord_equal() +
#      scale_fill_gradientn(colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
#                           na.value = "transparent",
#                           name = "°C x 10") +
#      labs(title = "Annual Mean Temperature",
#           x = "longitude",
#           y = "latitude") +
#      scale_x_continuous(expand = c(0, 0)) +
#      scale_y_continuous(expand = c(0, 0)) +
#      theme_minimal() +
#      theme(plot.title = element_text(hjust = 0.5),
#            axis.ticks.x = element_blank(),
#            axis.ticks.y = element_blank())

## ----load SDMtune-------------------------------------------------------------
#  library(SDMtune)

## ----help dataset-------------------------------------------------------------
#  help(virtualSp)
#  p_coords <- virtualSp$presence
#  bg_coords <- virtualSp$background

## ----plot presence------------------------------------------------------------
#  ggplot(map_data("world"), aes(long, lat)) +
#      geom_polygon(aes(group = group), fill = "grey95", color = "gray40", size = 0.2) +
#      geom_jitter(data = p_coords, aes(x = x, y = y), color = "red",
#                  alpha = 0.4, size = 1) +
#      labs(x = "longitude", y = "latitude") +
#      theme_minimal() +
#      theme(legend.position = "none") +
#      coord_fixed() +
#      scale_x_continuous(limits = c(-125, -32)) +
#      scale_y_continuous(limits = c(-56, 40))

## ----plot bg_model locations--------------------------------------------------
#  ggplot(map_data("world"), aes(long, lat)) +
#      geom_polygon(aes(group = group), fill = "grey95", color = "gray40", size = 0.2) +
#      geom_jitter(data = as.data.frame(bg_coords), aes(x = x, y = y),
#                  color = "blue", alpha = 0.4, size = 0.5) +
#      labs(x = "longitude", y = "latitude") +
#      theme_minimal() +
#      theme(legend.position = "none") +
#      coord_fixed() +
#      scale_x_continuous(limits = c(-125, -32)) +
#      scale_y_continuous(limits = c(-56, 40))

## ----prepare SWD--------------------------------------------------------------
#  data <- prepareSWD(species = "Virtual species", p = p_coords, a = bg_coords,
#                     env = predictors, categorical = "biome")
#  

## ----show SWD object----------------------------------------------------------
#  data

## ----show data----------------------------------------------------------------
#  head(data@data)

## ----show coords data---------------------------------------------------------
#  head(data@coords)

## ----show species-------------------------------------------------------------
#  data@species

## ----save SWD single, eval=FALSE----------------------------------------------
#  swd2csv(data, file_name = "data.csv")

## ----save SWDdouble, eval=FALSE-----------------------------------------------
#  swd2csv(data, file_name = c("presence.csv", "background.csv"))

## ----train--------------------------------------------------------------------
#  model <- train(method = "Maxent", data = data)

## ----print model--------------------------------------------------------------
#  model

## ----get slots----------------------------------------------------------------
#  slotNames(model)

## ----model slots--------------------------------------------------------------
#  slotNames(model@model)

## ----retrain------------------------------------------------------------------
#  model <- train(method = "Maxent", data = data, fc = "lqph", reg = 1, iter = 500)

## ----model witout default arguments-------------------------------------------
#  model <- train(method = "Maxent", data = data, fc = "lh", reg = 0.5, iter = 700)

## ----predict train------------------------------------------------------------
#  pred <- predict(model, data = data, type = "cloglog")

## ----print pred---------------------------------------------------------------
#  head(pred)

## ----predict presence---------------------------------------------------------
#  p <- data@data[data@pa == 1, ]
#  pred <- predict(model, data = p, type = "cloglog")
#  tail(pred)

## ----predict raster-----------------------------------------------------------
#  map <- predict(model, data = predictors, type = "cloglog")

## ----print raster output------------------------------------------------------
#  map

## ----save map, eval=FALSE-----------------------------------------------------
#  map <- predict(model, data = predictors, type = "cloglog", file = "my_map",
#                 format = "GTiff")

## ----exercise-----------------------------------------------------------------
#  # First create the extent that surrounds Cile
#  e = raster::extent(c(-77, -60, -56, -15))
#  # Now use the extent to make the prediction
#  map_e <- predict(model, data = predictors, type = "cloglog", extent = e)

## ----plot map default---------------------------------------------------------
#  plotPred(map)

## ----plot map custom----------------------------------------------------------
#  plotPred(map, lt = "Habitat\nsuitability",
#           colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

## ----thresholds---------------------------------------------------------------
#  ths <- thresholds(model, type = "cloglog")
#  ths

## ----plot pa------------------------------------------------------------------
#  plotPA(map, th = ths[3, 2])

## ----plot and save pa, eval=FALSE---------------------------------------------
#  plotPA(map, th = ths[3, 2], filename = "my_pa_map", format = "GTiff")

## ----auc----------------------------------------------------------------------
#  auc(model)

## ----plot roc-----------------------------------------------------------------
#  plotROC(model)

## ----tss----------------------------------------------------------------------
#  tss(model)

## ----aicc---------------------------------------------------------------------
#  aicc(model, env = predictors)

## ----train test---------------------------------------------------------------
#  library(zeallot)  # For unpacking assignment
#  c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE,
#                                   seed = 25)

## ----trai with train dataset--------------------------------------------------
#  model <- train("Maxent", data = train)

## ----evaluate test------------------------------------------------------------
#  auc(model)
#  auc(model, test = test)

## ----plot AUC train/test------------------------------------------------------
#  plotROC(model, test = test)

## ----experiment---------------------------------------------------------------
#  # Create an empty data.frame
#  output <- data.frame(matrix(NA, nrow = 10, ncol = 3))
#  colnames(output) <- c("seed", "trainAUC", "testAUC")
#  # Create 10 different random seeds
#  set.seed(25)
#  seeds <- sample.int(1000, 10)
#  
#  # Loop through the seeds
#  for (i in 1:length(seeds)) {
#    # Make the train/test split
#    c(train, test) %<-% trainValTest(data, test = 0.2, seed = seeds[i],
#                                     only_presence = TRUE)
#    # train the model
#    m <- train("Maxent", data = train)
#    # Populate the output data.frame
#    output[i, 1] <- seeds[i]
#    output[i, 2] <- auc(m)
#    output[i, 3] <- auc(m, test = test)
#  }
#  
#  # Print the output
#  output
#  # compute the range of the testing AUC
#  range(output[, 3])

## ----random folds-------------------------------------------------------------
#  folds <- randomFolds(data, k = 4, only_presence = TRUE, seed = 25)

## ----cv, eval=FALSE-----------------------------------------------------------
#  cv_model <- train("Maxent", data = data, folds = folds)
#  cv_model

## -----------------------------------------------------------------------------
#  auc(cv_model)
#  auc(cv_model, test = TRUE)

## ----enmeval block------------------------------------------------------------
#  library(ENMeval)
#  block_folds <- get.block(occ = data@coords[data@pa == 1, ],
#                           bg.coords = data@coords[data@pa == 0, ])
#  model <- train(method = "Maxent", data = data, fc = "l", reg = 0.8,
#                 folds = block_folds)

## ----enmeval checherboard1----------------------------------------------------
#  cb_folds <- get.checkerboard1(occ = data@coords[data@pa == 1, ],
#                                env = predictors,
#                                bg.coords = data@coords[data@pa == 0, ],
#                                aggregation.factor = 4)
#  model <- train(method = "Maxent", data = data, fc = "l", reg = 0.8,
#                 folds = cb_folds)

## ----blockCV------------------------------------------------------------------
#  library(blockCV)
#  library(raster)
#  # Create spatial points data frame
#  sp_df <- SpatialPointsDataFrame(data@coords,
#                                  data = as.data.frame(data@pa),
#                                  proj4string = crs(predictors))
#  e_folds <- envBlock(rasterLayer = predictors, speciesData = sp_df,
#                      species = "data@pa", k = 4, standardization = "standard",
#                      rasterBlock = FALSE, numLimit = 100)
#  model <- train(method = "Maxent", data = data, fc = "l", reg = 0.8,
#                 folds = e_folds)

