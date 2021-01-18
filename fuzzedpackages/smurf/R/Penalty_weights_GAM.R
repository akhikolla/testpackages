###############################################
#
# GAM coefficients used for penalty weights
#
###############################################

# formula: A formula object describing the model to be fitted
# family: A family object specifying the error distribution and link function for the model
# data: A data frame containing the model response and predictors for n observations
# weights: Vector of prior weights
# offset: Vector containing the offset for the model
# n.par.cov: List with number of parameters to estimate per predictor
# refcat: Indicator for reference category for Fused Lasso, Generalized Fused Lasso and Graph-Guided Fused Lasso
.gam.coef.pen.weights <- function(formula, family, data, weights, offset, n.par.cov, refcat) {
  
  # Phoney model frame without dropped levels to get penalty attributes (which can be dropped 
  # with drop.unused.levels = TRUE)
  mf_nodrop <- model.frame(formula = terms(formula, specials = "p", data = data), 
                           data = data, drop.unused.levels = FALSE)
  
  # Create new formula
  formula.gam.pw <- paste(as.character(formula[2]), "~ ")
  data.gam <- data
  
  # Add intercept to formula if present
  if (!attr(attr(mf_nodrop, "terms"), "intercept")) {
    # formula.gam.pw <- paste(formula.gam.pw, "-1")
    stop("GAM penalty weights are only computed if an intercept is included.")
    
  } else {
    formula.gam.pw <- paste(formula.gam.pw, "1")
  }
  
  for (j in 2:ncol(mf_nodrop)) {

    if (is.null(attr(mf_nodrop[, j], "cov.name"))) {
      # No covariate name given, use column name of mf_nodrop
      cov.name <- colnames(mf_nodrop)[j]
      refcat.cov <- ifelse(is.factor(mf_nodrop[, j]), 1L, 0L)
      
    } else {
      # Get covariate name
      cov.name <- attr(mf_nodrop[, j], "cov.name")
      refcat.cov <- attr(mf_nodrop[, j], "refcat")
    }
    
   

    # No penalty type given, set type to "none"
    if (is.null(attr(mf_nodrop[, j], "penalty"))) {
      pen.cov.gam <- "none"
      
    } else {
      # Get penalty type
      pen.cov.gam <- attr(mf_nodrop[, j], "penalty")
    }
    
    # Only relevel if no penalty, or if reference category present and one of the specified penalties 
    if (pen.cov.gam == "none" | (pen.cov.gam %in% c("flasso", "gflasso", "ggflasso") & refcat)) {
      
      if (is.factor(data.gam[, cov.name]) & refcat.cov > 1) {
        # Relevel to make chosen category the reference category
        data.gam[, cov.name] <- relevel(data.gam[, cov.name], ref = levels(data.gam[, cov.name])[refcat.cov])
      }
    }
    
    # Special handling of none, lasso, grouplasso and fused lasso
    if (pen.cov.gam %in% c("none", "lasso", "grouplasso", "flasso")) {
      
     
      if (is.factor(data.gam[, cov.name])) {
        # Categorical predictor
        
        if (any(is.na(suppressWarnings(as.numeric(levels(data.gam[, cov.name])))))) {
          # Non-numeric factor levels
          pen.cov.gam <- paste0(pen.cov.gam, ".factor.nonnumeric")
          
        } else {
          # Numeric factor levels
          pen.cov.gam <- paste0(pen.cov.gam, ".factor.numeric")
        }
        
      } else {
        # Lasso or Group Lasso with continuous predictor
        pen.cov.gam <- paste0(pen.cov.gam, ".numeric")
      }
    }
    
    if (pen.cov.gam %in% c("lasso.factor.nonnumeric", "grouplasso.factor.nonnumeric")) {
      stop("GAM penalty weights do not work when a non-numeric factor has a Lasso or Group Lasso penalty. Please use GLM penalty weights instead.")
    }
    
    if (pen.cov.gam == "gflasso" & !refcat) {
      stop("GAM penalty weights do not work for 'gflasso' and lambda1 > 0 or lambda2 > 0. Please use GLM penalty weights instead.")
    }
    
    if (pen.cov.gam == "flasso.factor.nonnumeric" & !refcat) {
      stop("GAM penalty weights do not work when a non-numeric factor has a Fused Lasso penalty and lambda1 > 0 or lambda2 > 0. Please use GLM penalty weights instead.")
    }
    
    # Add predictors to formula	
    if (pen.cov.gam %in% c("none.numeric", "none.factor.nonnumeric", 
                           "lasso.numeric", "lasso.factor.nonnumeric",
                           "grouplasso.numeric", "grouplasso.factor.nonnumeric", 
                           "flasso.factor.nonnumeric", "gflasso")) {
      
      # Add covariate to formula (no smooth effect)
      formula.gam.pw <- paste0(formula.gam.pw, " + ", cov.name)
      
    } else if (pen.cov.gam %in% c("none.factor.numeric", "lasso.factor.numeric", 
                                  "grouplasso.factor.numeric", "flasso.factor.numeric")) {
      
     
      if (pen.cov.gam %in% c("lasso.factor.numeric", "grouplasso.factor.numeric") | 
          (pen.cov.gam == "flasso.factor.numeric" & !refcat)) {
        # No reference class
        pc <- NULL
        
      } else {
        
        # Reference class (first level)
        pc <- levels(data.gam[, cov.name])[1]
      }

      # Add smooth effect to formula 
      formula.gam.pw <- paste0(formula.gam.pw, " + s(", paste0(cov.name, ".numer"), 
                               ", pc = ", pc, ") ")
      
      # Convert to numeric
      data.gam[, paste0(cov.name, ".numer")] <- as.numeric(levels(data.gam[, cov.name]))[data.gam[, cov.name]]
    
    } else if (pen.cov.gam == "2dflasso") {
      
      # Get covariate names
      cov.names <- attr(mf_nodrop[, j], "cov.names")
      cov.names <- c(cov.names[[1]][1], cov.names[[2]][1])
      # Remove ".binned" if present
      cov.names <- gsub(".binned*", "", cov.names)
      
      # Check if cov.names exist as predictor names
      if (any(!cov.names %in% colnames(data.gam))) {
        stop("Binned factors for interaction should have the original predictor name + '.binned' as predictor names, e.g. 'age.binned.")
      }
      
      # Create name of dummy
      dummy.name <- paste0("dummy.", cov.names[1], ".", cov.names[2])
      # Reference levels
      lev.ref <- as.numeric(c(levels(data.gam[, cov.names[1]])[1], levels(data.gam[, cov.names[2]])[1]))
      # Dummy for not being in the reference class
      dummy.gam <- rep(1, nrow(data.gam))
      dummy.gam[data.gam[, cov.names[1]] == lev.ref[1] | data.gam[, cov.names[2]] == lev.ref[2]] <- 0

      
      # Add to data
      data.gam[, dummy.name] <- dummy.gam 
      
      # Convert covariates to numeric and add to data frame if not already present
      if (!(paste0(cov.names[1], ".numer") %in% names(data.gam))) {
        l1 <- suppressWarnings(as.numeric(levels(data.gam[, cov.names[1]]))[data.gam[, cov.names[1]]])
        
        if (any(is.na(l1))) {
          stop("GAM penalty weights do not work when a non-numeric factor is included in the interaction. Please use GLM penalty weights instead.")
          
        } else {
          data.gam[, paste0(cov.names[1], ".numer")] <- l1
        }
      }
      
      if (!(paste0(cov.names[2], ".numer") %in% names(data.gam))) {
        l2 <- suppressWarnings(as.numeric(levels(data.gam[, cov.names[2]]))[data.gam[, cov.names[2]]])
        
        if (any(is.na(l2))) {
          stop("GAM penalty weights do not work when a non-numeric factor is included in the interaction. Please use GLM penalty weights instead.")
          
        } else {
          data.gam[, paste0(cov.names[2], ".numer")] <- l2
        }
      }
      
      # Add smooth interaction with dummy
      formula.gam.pw <- paste0(formula.gam.pw, " + s(", paste0(cov.names[1], ".numer"), ", ", paste0(cov.names[2], ".numer"),  
                               ", by = ", dummy.name, ")")
     
    } else if (pen.cov.gam == "ggflasso") {
      
      # Add .lat and .lon
      cov.names <- c(paste0(cov.name, ".lat"), paste0(cov.name, ".lon"))
      
      
      # Index of reference class, only used if refcat=TRUE
      ref.ind <- which(data.gam[, cov.name] == levels(data.gam[, cov.name])[1])[1]
      
      # Check if cov.names exist as predictor names
      if (!all(cov.names %in% colnames(data.gam))) {
        warning(paste0("Latitude and longitude values for 'ggflasso' should have the original predictor name", 
                        " + '.lat' and '.lon' as predictor names, e.g. 'zip.lat' and 'zip.lon. Since they are not given, ",
                        "a non-smooth effect is used for predictor ", deparse(substitute(cov.name)), "."))
        
        if (!refcat) {
          stop("GAM penalty weights do not work for 'ggflasso' and lambda1 > 0 or lambda2 > 0 when no latitude and longitude are given. Please use GLM penalty weights instead.")
        }
        # Add covariate to formula (no smooth effect)
        formula.gam.pw <- paste0(formula.gam.pw, " + ", cov.name)
        
      } else {
        
        if (refcat) {
          # Add smooth interaction with longitude and latitude
          formula.gam.pw <- paste0(formula.gam.pw, " + s(", cov.names[1], ", ", cov.names[2], 
                                   ", pc = c(", data.gam[ref.ind, cov.names[1]], ", ", data.gam[ref.ind, cov.names[2]], "))") 
        } else {
          # Add smooth interaction with longitude and latitude
          formula.gam.pw <- paste0(formula.gam.pw, " + s(", cov.names[1], ", ", cov.names[2], ")") 
        }

      }

    } else {
      stop("This penalty type is not handled correctly for GAM penalty weights.")
    }
  }
  
  # Fit GAM model
  gam.fit.pw <- gam(formula = as.formula(formula.gam.pw), data = data.gam, family = family, 
                    weights = weights, offset = offset)
  
  
  # Coefficient vector
  beta <- numeric(sum(sapply(n.par.cov, sum)))
  
  # Index for coefficient vector
  ind <- 0L
  
  pen.cov.gam <- vector("list", length(n.par.cov))
  
  if (attr(attr(mf_nodrop, "terms"), "intercept")) {
    
    # Predict GAM model to get intercept
    pg <- predict.gam(gam.fit.pw, newdata = data.gam[1,], type = "terms")
    
    beta[1] <- attr(pg, "constant")
    
    pen.cov.gam[[1]] <- "intercept"
    
    ind <- 1L
  }
  
  for (j in 2:ncol(mf_nodrop)) {
    
    if (is.null(attr(mf_nodrop[, j], "cov.name"))) {
      # No covariate name given, use column name of mf_nodrop
      cov.name <- colnames(mf_nodrop)[j]

    } else {
      # Get covariate name
      cov.name <- attr(mf_nodrop[, j], "cov.name")
    }

    
    # No penalty type given, set type to none
    if (is.null(attr(mf_nodrop[, j], "penalty"))) {
      pen.cov.gam[[j]] <- "none"
      
    } else {
      # Get penalty type
      pen.cov.gam[[j]] <- attr(mf_nodrop[, j], "penalty")
    }
    
    # Add predictors to formula
    if (pen.cov.gam[[j]] %in% c("none", "lasso", "grouplasso")) {
      
      if (is.factor(data.gam[, cov.name])) {
        # Categorical covariate
        
        if (any(is.na(suppressWarnings(as.numeric(levels(data.gam[, cov.name])))))) {
          # No numeric levels
          
          # Predict GAM model for current predictor
          pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", terms = cov.name)
          
        } else {
          # Numeric levels
          
          # Predict GAM model for current predictor
          pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", 
                            terms = paste0("s(", cov.name, ".numer)"))
        }
        
        if (pen.cov.gam[[j]] == "none") {
          # With reference level
          
          # Get coefficient per level (except reference level)
          for (i in 1:(length(levels(data.gam[, cov.name]))-1L)) {
            ind.none <- which(data.gam[, cov.name] == levels(data.gam[, cov.name])[i+1L])[1]
            # Avoid problems if level is not present
            if (length(ind.none) > 0) {
              beta[ind + i] <- pg[ind.none]
            }
          }
          
          ind <- ind + length(levels(data.gam[, cov.name])) - 1L
          
        } else {
          # No reference level
          
          # Get coefficient per level
          for (i in 1:length(levels(data.gam[, cov.name]))) {
            ind.lasso <- which(data.gam[, cov.name] == levels(data.gam[, cov.name])[i])[1]
            # Avoid problems if level is not present
            if (length(ind.lasso) > 0)
              beta[ind + i] <- pg[ind.lasso]
          }
          
          ind <- ind + length(levels(data.gam[, cov.name]))
        }
        
        
        
      } else {
        # Continuous numeric covariate

        # Get coefficient
        beta[ind + 1L] <- coef(gam.fit.pw)[cov.name]
        
        ind <- ind + 1L
      }

      
    } else if (pen.cov.gam[[j]] == "flasso") {
      
      if (any(is.na(suppressWarnings(as.numeric(levels(data.gam[, cov.name])))))) {
        # No numeric levels
        
        # Predict GAM model for current predictor
        pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", terms = cov.name)
        
      } else {
        # Numeric levels
        
        # Predict GAM model for current predictor
        pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", 
                          terms = paste0("s(", cov.name, ".numer)"))
      }
      
      # Get coefficient per level (except reference level)
      for (i in 1:(length(levels(data.gam[, cov.name]))-refcat)) {
        ind.flasso <- which(data.gam[, cov.name] == levels(data.gam[, cov.name])[i+refcat])[1]
        # Avoid problems if level is not present
        if (length(ind.flasso) > 0) {
          beta[ind + i] <- pg[ind.flasso]
        }
      }
      
      ind <- ind + length(levels(data.gam[, cov.name])) - refcat
      
    } else if (pen.cov.gam[[j]] == "gflasso") {
      
      # Predict GAM model for current predictor
      pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", terms = cov.name)
      
      # Get coefficient per level (except reference level)
      for (i in 1:(length(levels(data.gam[, cov.name]))-refcat)) {
        ind.gflasso <- which(data.gam[, cov.name] == levels(data.gam[, cov.name])[i+refcat])[1]
        # Avoid problems if a level is not present
        if (length(ind.gflasso) > 0) {
          beta[ind + i] <- pg[ind.gflasso]
        }
        
      }
      
      ind <- ind + length(levels(data.gam[, cov.name])) - refcat
      
    } else if (pen.cov.gam[[j]] == "2dflasso") {
      
      # Get covariate names
      cov.names <- attr(mf_nodrop[, j], "cov.names")
      cov.names_orig <- c(cov.names[[1]][1], cov.names[[2]][1])
      # Remove ".binned" if present
      cov.names <- gsub(".binned*", "", cov.names_orig)
      
      # Predict GAM model for current predictor
      pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", 
                        terms = paste0("s(", cov.names[1], ".numer,", cov.names[2], ".numer):dummy.",
                                       cov.names[1], ".", cov.names[2]))

      # Levels of first (binned) predictor
      lev.1 <- levels(data.gam[, cov.names_orig[1]])
      # Levels of second (binned) predictor
      lev.2 <- levels(data.gam[, cov.names_orig[2]])
      
      ind2 <- 1L
      for (i in 1:(length(lev.1)-1L)) {
        for (k in 1:(length(lev.2)-1L)) {
          beta[ind + ind2] <- median(unique(pg[data.gam[, cov.names_orig[1]] == lev.1[i+1L] & 
                                            data.gam[, cov.names_orig[2]] == lev.2[k+1L]]))
          ind2 <- ind2 + 1L
        }
      }
      
      ind <- ind + (length(lev.1) - 1L) * (length(lev.2) - 1L)
      
    } else if (pen.cov.gam[[j]] == "ggflasso") {
    
      
      # Check if exist as predictor names
      if (!all(c(paste0(cov.name, ".lat"), paste0(cov.name, ".lon")) %in% colnames(data.gam))) {
        
        # Predict GAM model for current predictor
        pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", terms = cov.name)
      
      } else {
        
        # Predict GAM model for current predictor
        pg <- predict.gam(gam.fit.pw, newdata = data.gam, type = "terms", 
                          terms = paste0("s(", cov.name, ".lat,", cov.name, ".lon)"))
        
      }
      
      # Get coefficient per level (except reference level)
      for (i in 1:(length(levels(data.gam[, cov.name]))-refcat)) {
        ind.ggflasso <- which(data.gam[, cov.name] == levels(data.gam[, cov.name])[i+refcat])[1]
        # Avoid problems if a level is not present
        if (length(ind.ggflasso) > 0) {
          beta[ind + i] <- pg[ind.ggflasso, ]
        }
      }
      ind <- ind + length(levels(data.gam[, cov.name])) - refcat
    }
    
  }

  return(beta)
}