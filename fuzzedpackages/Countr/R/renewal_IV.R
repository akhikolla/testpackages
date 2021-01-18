## check initial values
.checkInitialValues <- function(dist, start, modelMatrixList, weights, Y, ...) {
    switch(dist,
           "weibull" = {
               Xscale <- modelMatrixList[["scale"]]
               nm <- gsub('\\(Intercept\\)', "", paste0("scale_", colnames(Xscale)))

               if ("shape" %in% names(modelMatrixList)) {
                   Xshape <- modelMatrixList[["shape"]]
                   nm_shape <- gsub('\\(Intercept\\)', "", paste0("shape_",
                                                                  colnames(Xshape)))
                   nm <- c(nm, nm_shape)
               } else
                   nm <- c(nm, "shape_")

               valid <- FALSE
               if (!is.null(start)) {
                   if(all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unnamed initial values found!",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!"shape" %in% names(modelMatrixList)) {
                       IV <- glm.fit(Xscale, Y, family = poisson(),
                                     weights = weights)
                       start <- c(IV$coefficients, log(1))
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters!"))
               }
           },
           "weibullgam" = {
               Xscale <- modelMatrixList[["scale"]]
               nm <- gsub('\\(Intercept\\)', "", paste0("scale_", colnames(Xscale)))
               if ("shape_" %in% names(modelMatrixList)) {
                   Xshape <- modelMatrixList[["shape"]]
                   nm_shape <- gsub('\\(Intercept\\)', "", paste0("shape_", colnames(Xshape)))
                   nm <- c(nm, nm_shape)
               } else
                   nm <- c(nm, "shape_")

               if ("shapeGam_" %in% names(modelMatrixList)) {
                   XshapeGam <- modelMatrixList[["shapeGam"]]
                   nm_shapeGam <- gsub('\\(Intercept\\)', "",
                                    paste0("shapeGam_", colnames(XshapeGam))
                                    )
                   nm <- c(nm, nm_shapeGam)
               } else
                   nm <- c(nm, "shapeGam_")

               if ("scaleGam_" %in% names(modelMatrixList)) {
                   XscaleGam <- modelMatrixList[["scaleGam_"]]
                   nm_sclaeGam <- gsub('\\(Intercept\\)', "",
                                       paste0("scaleGam_", colnames(XscaleGam)))
                   nm <- c(nm, nm_sclaeGam)
               } else
                   nm <- c(nm, "scaleGam_")

               valid <- FALSE
               if (!is.null(start)) {
                   if(all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unnamed initial values found!",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (length(modelMatrixList) == 1) {
                       IV <- glm.fit(Xscale, Y, family = poisson(), weights = weights)
                       ## change parametrization in terms of 1/r and 1/alpha (the gamma pars)
                       start <- c(IV$coefficients, 1, 16, 4)
                       names(start) <- nm
                   } else
                      stop(paste("initial values should be provided in control",
                                 "when regression on ancillary parameters!")
                           )
               }
           },
           "gamma" = {
               Xrate <- modelMatrixList[["rate"]]
               nm <- gsub('\\(Intercept\\)', "",
                          paste0("rate_", colnames(Xrate)))

               if ("shape" %in% names(modelMatrixList)) {
                   Xshape <- modelMatrixList[["shape"]]
                   nm_shape <- gsub('\\(Intercept\\)', "",
                                    paste0("shape_", colnames(Xshape)))
                   nm <- c(nm, nm_shape)
               } else
                   nm <- c(nm, "shape_")

               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unnamed initial values found!",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!"shape" %in% names(modelMatrixList)) {
                       IV <- glm.fit(Xrate, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, log(1))
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters!"))
               }
           },
           "gengamma" = {
               Xmu <- modelMatrixList[["mu"]]
               nm <- paste0("mu_", colnames(Xmu))

               if ("sigma" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("sigma_",
                                      colnames(modelMatrixList[["sigma"]])))
               else
                   nm <- c(nm, "sigma_")

               if ("Q" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("Q_",
                                      colnames(modelMatrixList[["Q"]])))
               else
                   nm <- c(nm, "Q_")

               nm <- gsub('\\(Intercept\\)', "", nm)
               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unnamed initial values found!",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!any(c("sigma" %in%  names(modelMatrixList),
                              "Q" %in%  names(modelMatrixList))
                            )
                       ) {
                       IV <- glm.fit(Xmu, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, log(1), 1)
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters!"))
               }
           },
           "burr" = {
               X <- modelMatrixList[["scale"]]
               nm <- paste0("scale_", colnames(X))

               if ("shape1" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("shape1_",
                                      colnames(modelMatrixList[["shape1"]])))
               else
                   nm <- c(nm, "shape1_")

               if ("shape2" %in% names(modelMatrixList))
                   nm <- c(nm, paste0("shape2_",
                                      colnames(modelMatrixList[["shape2"]])))
               else
                   nm <- c(nm, "shape2_")

               nm <- gsub('\\(Intercept\\)', "", nm)
               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unnamed initial values found!",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid) {
                   if (!any(c("shape1" %in%  names(modelMatrixList),
                              "shape2" %in%  names(modelMatrixList))
                            )
                       ) {
                       IV <- glm.fit(X, Y, family = poisson(), weights = weights)
                       start <- c(IV$coefficients, log(1), log(3))
                       names(start) <- nm
                   } else
                       stop(paste("initial values should be provided in control",
                                  "when regression on ancillary parameters!"))
               }
           },
           "custom" = {
               customPars <- list(...)[[1]]
               parNames <- customPars[["parNames"]]

               location <- parNames[1]
               X <- modelMatrixList[[location]]
               nm <- paste0(location, "_", colnames(X))

                if (length(parNames) > 1) {
                   for (i in 2:length(parNames)) {
                       pari <- parNames[i]
                       if (pari %in% names(modelMatrixList)) {
                           nm <- c(nm, paste0(pari, "_",
                                              colnames(modelMatrixList[[pari]])
                                              )
                                   )
                       } else
                           nm <- c(nm, paste0(pari, "_"))
                   }
                }

               nm <- gsub('\\(Intercept\\)', "", nm)
               valid <- FALSE
               if (!is.null(start)) {
                   if (all(names(start) == nm))
                       valid <- TRUE
                   else if (length(start) == length(nm) &
                            all(!is.na(start))) {
                       warning(paste("unnamed initial values found!",
                                     "They will be used as they are")
                               )
                       names(start) <- nm
                       valid <- TRUE
                   }
               }

               if (!valid)
                   stop(paste("initial values should be provided in control",
                              "when custom survival functions are passed!"))
           }
           )
    start
}
