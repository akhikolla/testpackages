    # --------------------------------------------------------------------------------------------------------------------------
    ############################################## Estimations and tests functions ##############################################
    ## --------------------------------------------------------------------------------------------------------------------------

  # Estimation function of the mean number of mutations. Compute also the probability of mutation if with.prob is TRUE and the fitness parameter if fitness is empty.
    # Arguments:
    # mc:  a (non-empty) numeric vector of mutants counts.
    # fn:  an optional (non-empty) numeric vector with same length as mc of final numbers of cells.
    # mfn:  mean final number of cells which. Ignored if fn is non-missing.
    # cvfn:  coefficient of variation of final number of cells. Ignored if fn is non-missing.
    # fitness:  fitness parameter: ratio of growth rates of normal and mutant cells. If fitness is NULL (default) then the fitness will be estimated. Otherwise, the given value will be used to estimate the mean mutation number mutations
    # death:  death probability. Must be smaller than 0.5.
    # plateff: plating efficiency.
    # method:  estimation method as a character string: one of ML (default), P0, or GF.
    # winsor:  winsorization parameter: positive integer. Only used when method is ML or when method is P0 and fitness is NULL.
    # model:  statistical lifetime model as a character string: one of LD (default) for Luria-Delbruck model with exponential lifetimes, or H for Haldane model with constant lifetimes.
    # Assumptions for the estimation
    #   - plating efficiency
    #   - fitness (if given)
    #   - death probability
    #   - lifetimes distribution model: exponentialy distributed lifetime ("LD") or constant lifetime "H"
    #   - sample of mutants counts
    #   - sample of finals counts (or a mean and coefficient of variation for finals counts)
    # Three possible methods:
    #   - p0-method ("P0")   (winsorize the mc sample if rho is estimated with the parameter winsor)
    #   - Generating function method ("GF")
    #   - Maximum of Likelihood method ("ML")  (winsorize the mc sample with the parameter winsor)
    #
    # Returns:
    #   - Estimate of alpha (or pi) and rho if non given
    #   - Standard deviation of alpha, (or pi) and rho if non given

mutestim <- function(mc, fn = NULL, mfn = NULL, cvfn = NULL,                 # user's data
                  fitness = NULL, death = 0., plateff = 1.,                # user's parameters
                  model = c("LD", "H", "I"), muinf = +Inf,                       # clone growth model
                  method = c("ML", "GF", "P0"), winsor = 2000) # estimation method
                  {

    if(missing(method)){method <- "ML"}
    if(missing(model)){model <- "LD"}
    method <- match.arg(method, several.ok = TRUE)
    if(length(method) > 1){
      if(sum(method == method[1]) == length(method)) method <- method[1]
      else stop("If you use a 'data.table', 'method' can not have different values in a same class.")
    }
    model <- match.arg(model, several.ok = TRUE)
    if(length(model) > 1){
      if(sum(model == model[1]) == length(model)) model <- model[1]
      else stop("If you use a 'data.table', 'model' can not have different values in a same class.")
    }

    if(is.null(mc)){
      stop("'mc' is empty...")
    }
    # if(!is.null(fn)){
    #   if(!is.null(mfn) | !is.null(cvfn)){
    #     warning("When 'fn' is non-empty, 'mfn' and 'cvfn' are ignored.")
    #   }
    # }
    if(sum(mc == 0) == length(mc)){
      warning("'mc' does not contain mutants count...")
      output <- list(mutations = 0., sd.mutations = 0., fitness = 1., sd.fitness = 0)
      return(output)
    }
    if(sum(floor(mc) == mc) != length(mc)) stop("'mc' must be a vector of integer.")
    if(!is.null(fn)){
      if(length(fn) != length(mc)){
	stop("'fn must have the same length as 'mc'.")
      }
      mfn <- mean(fn)
      cvfn <- sd(fn)/mfn
      if(sd(fn) == 0) fn <- NULL
    }
    if(!is.null(mfn)){
      if(length(mfn) > 1){
        if(sum(mfn == mfn[1]) == length(mfn)) mfn <- mfn[1]
        else stop("If you use a 'data.table', 'mfn' can not have different values in a same class.")
      }
      if(mfn < 0){
	stop("'mfn' must be a positive number.")
      }
      if(is.null(cvfn)) cvfn <- 0
    }
    if(!is.null(cvfn)){
      if(length(cvfn) > 1){
        if(sum(cvfn == cvfn[1]) == length(cvfn)) cvfn <- cvfn[1]
        else stop("If you use a 'data.table', 'cvfn' can not have different values in a same class.")
      }
      if(cvfn < 0){
	stop("'cvfn' must be a positive number.")
      }
      if(is.null(mfn)) mfn <- 1e9
    }
    if(length(fitness) > 1){
      if(sum(fitness == fitness[1]) == length(fitness)) fitness <- fitness[1]
      else stop("If you use a 'data.table', 'fitness' can not have different values in a same class.")
    }
    if(!is.null(fitness)){
      if(fitness < 0){
	stop("'fitness' must be empty or a positive number.")
      }
    }
    if(length(death) > 1){
      if(sum(death == death[1]) == length(death)) death <- death[1]
      else stop("If you use a 'data.table', 'death' can not have different values in a same class.")
    }
    if(death < 0 | death >= 0.5){
      stop("'death must be a positive and < 0.5 number.")
    }
    if(length(plateff) > 1){
      if(sum(plateff == plateff[1]) == length(plateff)) plateff <- plateff[1]
      else stop("If you use a 'data.table', 'plateff' can not have different values in a same class.")
    }
    if(plateff > 1){
      stop("'plateff' must be a positive and <= 1 number.")
    }
    if(length(winsor) > 1){
      if(sum(winsor == winsor[1]) == length(winsor)) winsor <- winsor[1]
      else stop("If you use a 'data.table', 'winsor' can not have different values in a same class.")
    }
    if(winsor < 0 | trunc(winsor) != winsor){
      stop("'winsor' must be a single positive integer.")
    }
    if(model == "I"){
      if(!is.finite(muinf)){
      	warning("'muinf' is infinite, 'model' is set to 'LD'.")
      	model <- "LD"
      }
    } else {
      if(is.finite(muinf)) {
      	warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      	muinf <- +Inf
      }
    }
    if(method == "P0" & sum(mc == 0) == 0 & death == 0){
      stop("If 'death' is zero, P0 method can not be used if 'mc' does not contain any null counts.")
    }
    if(method == "P0" & plateff < 1){
      warning("'plateff' can not be taking into account for P0 method: 'plateff' is set to 1.")
      plateff <- 1
    }
    if(method == "ML" & model == "H" & plateff < 1){
      warning("'plateff' can not be taking into account for ML method when 'model' is 'H': 'plateff' is set to 1.")
      plateff <- 1
    }
    if(is.null(fitness)){   # If fitness is empty: compiute estimates of mean number of mutations (mutation probaility), fitness

      # Maximum of Likelihood estimators
      if(method == "ML"){
      	if(!is.null(fn)) output <- MutationProbabilityFitnessMLOptimization(mc = mc, fn = fn, death = death, plateff = plateff, model = model, winsor = winsor, muinf = muinf)
      	else output <- MutationFitnessMLOptimization(mc = mc, mfn = mfn, cvfn = cvfn, death = death, plateff = plateff, model = model, winsor = winsor, muinf = muinf)
      }
      # P0 method
      if(method == "P0") output <- MutationFitnessP0Estimation(mc = mc, fn = fn, mfn = mfn, cvfn = cvfn, death = death, model = model, winsor = winsor)

      # GF method
      if(method == "GF") {
      	output <- MutationFitnessGFEstimation(mc = mc, mfn = mfn, cvfn = cvfn, death = death, plateff = plateff, model = model, muinf = muinf, init = FALSE)
      	if(!output$succeeds) warning(paste("Impossible to estimate 'fitness' with 'GF'-method: 'fitness' is set to default value 1 and only", if(!is.null(mfn)){"mutation probability"}else{"mutation number"}, "is estimated.", sep = " "))
      	output$succeeds <- NULL
      }
    } else {    # Else: compute estimate(s) of mean number of mutations, or mutation probability

      # Maximum of Likelihood estimator of mutations or mutprob
      if(method == "ML"){
      	if(!is.null(fn)) output <- MutationProbabilityMLOptimization(mc = mc, fn = fn, fitness = fitness, death = death, plateff = plateff, model = model, winsor = winsor, muinf = muinf)
      	else output <- MutationMLOptimization(mc = mc, mfn = mfn, cvfn = cvfn, fitness = fitness, death = death, plateff = plateff, model = model, winsor = winsor, muinf = muinf)
      }
      # P0 method
      if(method == "P0") output <- MutationP0Estimation(mc, fn = fn, mfn = mfn, cvfn = cvfn, death = death)     # P0 estimator of mutations

      # GF method
      if(method == "GF")  output <- MutationGFEstimation(mc = mc, mfn = mfn, cvfn = cvfn, fitness = fitness, death = death, plateff = plateff, model = model, muinf = muinf, init = FALSE)

    }
    # if(plateff < 1) output <- lapply(output, function(o) o*(plateff-1)/(plateff*log(plateff)))
    output
}


    # One-sample or two-sample tests for mean number of mutations, fitness and mutation probability using mutestim estimates
    # The possible assumptions are the same as for the mutestim function.
    # Arguments:
    # Returns a "flantest" object, which contains:
    # Tstat:  the value of the computed statistic.
    # parameter:  the values of the parameter of the model: fitness(if not tested), death, mfn (if needed) and cvfn
    # p.value:  the p-value(s) of the test.
    # conf.int:  a confidence interval for the parameter(s) of interest appropriate to the specified alternative hypothesis.
    # estimates:  the estimate(s) of interest.
    # null.value:  the specified hypothesized value(s).
    # alternative:  a character string describing the alternative hypothesis.
    # model:  the statistical lifetime model.
    # method:  method used to compute the estimation(s) of the parameter(s) of interest.
    # data.name:  a character string giving the name of the complete data.


flan.test <- function(mc, fn = NULL, mfn = NULL, cvfn = NULL,                      # user's data
               fitness = NULL, death = 0., plateff = 1.,                        # user's parameters
               model = c("LD", "H", "I"), muinf = +Inf,                            # clone growth model
               mutations0 = 1., mutprob0 = NULL, fitness0 = 1.,       # null hypotheses
               conf.level = 0.95,                              # confidence level
               alternative = c("two.sided", "less", "greater"),  # alternative
               method = c("ML", "GF", "P0"), winsor = 2000)          # estimation method
               {

   with.prob <- FALSE               # Boolean: if TRUE (if fn, mfn, or cvfn are given), mutprob is tested instead of mutations

   if(missing(method)) method <- "ML"
   if(missing(model)) model <- "LD"


   method <- match.arg(method)
   model <- match.arg(model)


   if(is.null(mc)){
     stop("'mc' is empty...")
   } else {
    if(is.list(mc)){
    if(length(mc) == 1){
      nsamples <- 1
      dname <- list(deparse(substitute(mc)))
      if(is.null(mc[[1]])) stop("'mc[[1]]' is empty...")
      if(is.list(fn)){
        if(length(fn) != nsamples) stop("'mc' and 'fn' must have the same length.")
        if(!is.null(fn[[1]])){
          dname <- c(dname, deparse(substitute(fn)))
          with.prob <- TRUE
        }
      } else {
          if(!is.null(fn)) stop("'fn' must have the same type as 'mc'.")
          fn <- list(fn)
        }
      }
      if(length(mc) == 2){
        nsamples <- 2
        dname <- list(c(paste(deparse(substitute(mc)), "1", sep = ""), paste(deparse(substitute(mc)), "2", sep = "")))
        if(is.null(mc[[1]])) stop("'mc[[1]]' is empty...")
        if(is.null(mc[[2]])) stop("'mc[[2]]' is empty...")

        if(is.list(fn)){
          if(length(fn) != nsamples) stop("'fn' must have the same length as 'mc'.")

          if(!is.null(fn[[1]])) {
            with.prob <- TRUE
            if(length(fn[[1]]) != length(mc[[1]])) stop("'fn[[1]]' must have the same length as 'mc[[1]]'.")
            # if(!is.null(mfn[1]) | !is.null(cvfn[1])) {
            #   if(abs(mean(fn[[1]])-mfn[1]) > 1e-6 & abs(sd(fn[[1]])/mean(fn[[1]])-cvfn[1]) > 1e-6) {
            #     warning("'fn[[1]]' is non-empty: 'mfn[1]' and 'cvfn[1]' will be ignored.")
            #     mfn[1] <- mean(fn[[1]])
            #     cvfn[1] <- sd(fn[[1]])/mfn[1]
            #   }
            # }
          #   if(is.null(fn[[2]])){
          #     if(is.null(mfn[2]) & is.null(cvfn[2])){
          #       warning("'fn[[2]]' is empty: empirical informations of 'fn[[1]]' are considered for sample 2.")
          #       if(is.null(mfn[2])) mfn <- list(rep(mean(fn[[1]]), 2))
          #       if(is.null(cvfn[2])) cvfn <- list(rep(sd(fn[[1]])/mfn[[1]]))
          #     }
          #   # fn <- list(NULL, NULL)
          # }
           if(!is.null(fn[[2]])) {
              if(length(fn[[2]]) != length(mc[[2]])) stop("'fn[[2]]' must have the same length as 'mc[[2]]'.")
              # if(!is.null(mfn[2]) | !is.null(cvfn[2])) {
              #   if(abs(mean(fn[[2]])-mfn[2]) > 1e-6 & abs(sd(fn[[2]])/mean(fn[[2]])-cvfn[2]) > 1e-6) {
              #     warning("'fn[[2]]' is non-empty: 'mfn[2]' and 'cvfn[2]' will be ignored.")
              #     mfn[2] <- mean(fn[[2]])
              #     cvfn[2] <- sd(fn[[2]])/mfn[2]
              #   }
              # }
              dname <- c(dname, list(c(paste(deparse(substitute(fn)), "1", sep = ""), paste(deparse(substitute(fn)), "2", sep = ""))))
              with.prob <- TRUE
              # mfn <- unlist(lapply(fn, mean))
              # cvfn <- unlist(lapply(fn, sd))/mfn
            }
          # } else {
          #   # if(!is.null(fn[[2]])) {
          #   #   if(is.null(mfn[1]) & is.null(cvfn[1])) {
          #   #     if(is.null(mfn[1])) mfn <- rep(mean(fn[[2]]), 2)
          #   #     if(is.null(cvfn[1])) cvfn <- sd(fn[[2]])/mfn
          #   #     warning("'fn[[1]]' is empty: empirical informations of 'fn[[2]]' are considered for sample 1.")
          #   #   }
          #   #   # if(!is.null(mfn[2]) | !is.null(cvfn[2])) {
          #   #   #   if(abs(mean(fn[[2]])-mfn[2]) > 1e-6 & abs(sd(fn[[2]])/mean(fn[[2]])-cvfn[2]) > 1e-6) {
          #   #   #       warning("'fn[[2]]' is non-empty: 'mfn[2]' and 'cvfn[2]' will be ignored.")
          #   #   #       mfn[2] <- mean(fn[[2]])
          #   #   #       cvfn[2] <- sd(fn[[2]])/mfn[2]
          #   #   #   }
          #   #   # }
          #   # }
          } else {
            if(!is.null(fn[[2]])) with.prob <- TRUE
          }
        } else {
          if(!is.null(fn)) stop("'fn' must have the same type as 'mc'.")
          fn <- list(fn, fn)
        }
        }
      } else {
      nsamples <- 1
      dname <- list(deparse(substitute(mc)))
      if(!is.null(fn)){
        if(is.list(fn)) stop("'fn' must have the same type as 'mc'.")
        if(length(fn) != length(mc)) stop("'fn' must have the same length as 'mc'.")
        dname <- c(dname, deparse(substitute(fn)))
        with.prob <- TRUE
        # if(!is.null(mfn) | !is.null(cvfn)){
        #   if(abs(mean(fn)-mfn) > 1e-6 & abs(sd(fn)/mean(fn)-cvfn) > 1e-6) {
        #     warning("'fn' is non-empty: 'mfn' and 'cvfn' will be ignored.")
        #     mfn <- mean(fn)
        #     cvfn <- sd(fn)/mfn
        #   }
        # }
      }
      mc <- list(mc)
      fn <- list(fn)
    }
    if(nsamples == 1){
      if(!is.null(mfn)){
        if(!is.null(fn[[1]])){
          warning("'fn' is non-empty: 'mfn' is ignored.")
          mfn <- mean(fn[[1]])
        } else {
          if(length(mfn) > 1) {
            if(is.list(mfn)) mfn <- unlist(mfn)
            if(sum(mfn == mfn[1]) == length(mfn)) mfn <- mfn[1]
            else stop("If you use a 'data.table', 'mfn' can not have different values in a same class.")
          }
          if(mfn < 0) stop("'mfn' must be empty or a positive number.")
          with.prob <- TRUE
        }
        if(is.null(cvfn)) {
          warning("'cvfn' is empty but 'mfn' is not: 'cvfn' is set to 0.")
          cvfn <- if(!is.null(fn[[1]])) sd(fn[[1]])/mfn else 0
        }
      } else {
        if(!is.null(fn[[1]])) mfn <- mean(fn[[1]])
      }
      if(!is.null(cvfn)){
        if(!is.null(fn[[1]])){
          warning("'fn' is non-empty: 'cvfn' is ignored.")
          cvfn <- sd(fn[[1]])/mfn
        } else {
          if(length(cvfn) > 1) {
            if(is.list(cvfn)) cvfn <- unlist(cvfn)
            if(sum(cvfn == cvfn[1]) == length(cvfn)) cvfn <- cvfn[1]
            else stop("If you use a 'data.table', 'cvfn' can not have different values in a same class.")
            }
            if(cvfn < 0) stop("'cvfn' must be empty or a positive number.")
            if(is.null(mfn)) mfn <- if(!is.null(fn[[1]])) mean(fn[[1]]) else 1e9
            with.prob <- TRUE
          }
        } else {
          if(!is.null(fn[[1]])) cvfn <- sd(fn[[1]])/mfn
        }
        if(!is.null(fitness)){
          if(is.list(fitness)) fitness <- unlist(fitness)
          if(length(fitness) > 1){
            if(sum(fitness == fitness[1]) == length(fitness)) fitness <- fitness[1]
            else stop("If you use a 'data.table', 'fitness' can not have different values in a same class.")
          }
          if(fitness < 0) stop("'fitness' must be empty or a positive number.")
          fitness0 <- NULL
        }
        if(with.prob & missing(mutprob0)) mutprob0 <- 1/mfn

        if(is.list(death)) death <- unlist(death)
        if(length(death) > 1){
          if(sum(death == death[1]) == length(death)) death <- death[1]
          else stop("If you use a 'data.table', 'death' can not have different values in a same class.")
        }
        if(death >= 0.5) stop("'death' must be a positive and < 0.5 number.")

        if(is.list(plateff)) plateff <- unlist(plateff)
        if(length(plateff) > 1){
          if(sum(plateff == plateff[1]) == length(plateff)) plateff <- plateff[1]
          else stop("If you use a 'data.table', 'plateff' can not have different values in a same class.")
        }
        if(plateff > 1) stop("'plateff' must be a positive and <= 1 number.")
        if(plateff < 1 & method == "P0"){
          warning("'plateff' can not be taking into account for P0 method: 'plateff' is set to 1.")
          plateff <- 1
        }
        if(model == "H" & plateff < 1 & method == "ML"){
          warning("'plateff' can not be taking into account for ML method when 'model' is 'H': 'plateff' is set to 1.")
          plateff <- 1
          }
    } else {
      if(!is.null(mfn)){
        if(!is.list(mfn)){
          if(sum(mfn < 0) != 0 | length(mfn) > 2) stop("if given, 'mfn' must be a vector or a list with size <= 2 of positive numbers.")
          if(length(mfn == 1)) mfn <- c(mfn, mfn)
          mfn <- list(mfn[1], mfn[2])
        }
        if(sum(unlist(mfn) < 0) != 0 | length(mfn) > 2) stop("if given, 'mfn' must be a vector or a list with size <= 2 of positive numbers.")
        if(!is.null(mfn[[1]]) & !is.null(fn[[1]])){
          warning("'fn[[1]]' is non-empty: 'mfn[[1]]' is ignored.")
          mfn[[1]] <- mean(fn[[1]])
          # cvfn[[1]] <- sd(fn[[1]])/mfn[[1]]
        }
        if(!is.null(mfn[[2]]) & !is.null(fn[[2]])){
          warning("'fn[[2]]' is non-empty: 'mfn[[2]]' is ignored.")
          mfn[[2]] <- mean(fn[[2]])
          # cvfn[[2]] <- sd(fn[[2]])/mfn[[2]]
        }
        if(is.null(mfn[[1]])) {
          if(is.null(fn[[1]])){
            if(!is.null(fn[[2]])){
              warning("'fn[[1]]' is empty: 'mfn[[1]]' is set as mean of 'fn[[2]]'.")
              mfn[[1]] <- mean(fn[[2]])
            }
          } else {
            mfn[[1]] <- mean(fn[[1]])
          }
        } else {
          if(is.null(cvfn[[1]])) {
            warning("'cvfn[[1]]' is empty but 'mfn[[1]]' is not: 'cvfn[[1]]' is set to 0.")
            cvfn[[1]] <- 0
          }
        }
        if(is.null(mfn[[2]])) {
          if(is.null(fn[[2]])){
            if(!is.null(fn[[1]])){
              warning("'fn[[2]]' is empty: 'mfn[[2]]' is set as mean of 'fn[[1]]'.")
              mfn[[2]] <- mean(fn[[1]])
            }
          } else {
            mfn[[2]] <- mean(fn[[2]])
          }
        } else {
          if(is.null(cvfn[[2]])) {
            warning("'cvfn[[2]]' is empty but 'mfn[[2]]' is not: 'cvfn[[2]]' is set to 0.")
            cvfn[[2]] <- 0
          }
        }
        with.prob <- TRUE
      } else {
        if(!is.null(fn[[1]])) mfn <- mean(fn[[1]])
        if(!is.null(fn[[2]])) mfn <- c(mfn, mean(fn[[2]]))
      }

      if(!is.null(cvfn)){
        if(!is.list(cvfn)){
          if(sum(cvfn < 0) != 0 | length(cvfn) > 2) stop("if given, 'cvfn' must be a vector or a list with size <= 2 of positive numbers.")
          if(length(cvfn) == 1) cvfn <- c(cvfn, cvfn)
          cvfn <- list(cvfn[1], cvfn[2])
        }
        if(sum(unlist(cvfn) < 0) != 0 | length(cvfn) > 2) stop("if given, 'cvfn' must be a vector or a list with size <= 2 of positive numbers.")
        if(!is.null(cvfn[[1]]) & !is.null(fn[[1]])){
          warning("'fn[[1]]' is non-empty: 'cvfn[[1]]' is ignored.")
          cvfn[[1]] <- sd(fn[[1]])/mean(fn[[1]])
        }
        if(!is.null(cvfn[[2]]) & !is.null(fn[[2]])){
          warning("'fn[[2]]' is non-empty: 'cvfn[[2]]' is ignored.")
          cvfn[[2]] <- sd(fn[[2]])/mean(fn[[2]])
        }
        if(is.null(cvfn[[1]])) {
          if(is.null(fn[[1]])){
            if(!is.null(fn[[2]])){
              warning("'fn[[1]]' is empty: 'cvfn[[1]]' is set as coef. of variation of 'fn[[2]]'.")
              cvfn[[1]] <- sd(fn[[2]])/mean(fn[[2]])
            }
          } else {
            cvfn[[1]] <- sd(fn[[1]])/mean(fn[[1]])
          }
        }
        if(is.null(cvfn[[2]])) {
          if(is.null(fn[[2]])){
            if(!is.null(fn[[1]])){
              warning("'fn[[2]]' is empty: 'cvfn[[2]]' is set as coef. of variation of 'fn[[1]]'.")
              cvfn[[2]] <- sd(fn[[1]])/mean(fn[[1]])
            }
          } else {
            cvfn[[2]] <- sd(fn[[2]])/mean(fn[[2]])
          }
        }
        with.prob <- TRUE
      } else {
        if(!is.null(fn[[1]])) cvfn <- sd(fn[[1]])/mean(fn[[1]])
        if(!is.null(fn[[2]])) cvfn <- c(cvfn, sd(fn[[2]])/mean(fn[[2]]))
      }

      if(!is.null(fn[[1]])){
        if(is.null(fn[[2]])){
          if(is.null(mfn[2])) {
            warning("'fn[[2]]' is empty: 'mfn[[2]]' is set as mean of 'fn[[1]]'.")
            # mfn <- list(NULL, mean(fn[[1]]))
            mfn <- rep(mean(fn[[1]]), 2)
          }
          if(is.null(cvfn[2])) {
            warning("'fn[[2]]' is empty: 'cvfn[[2]]' is set as coef. of variation of 'fn[[1]]'.")
            # cvfn <- list(NULL, sd(fn[[1]])/mfn[[2]])
            cvfn <- rep(sd(fn[[1]])/mean(fn[[1]]), 2)
          }
        }
      } else {
        if(!is.null(fn[[2]])){
          if(is.null(mfn[[1]])){
            warning("'fn[[1]]' is empty: 'mfn[[1]]' is set as mean of 'fn[[2]]'.")
            # mfn <- list(mean(fn[[2]]), NULL)
            mfn <- rep(mean(fn[[2]]), 2)
          }
          if(is.null(cvfn[1])) {
            warning("'fn[[1]]' is empty: 'cvfn[[1]]' is set as coef. of variation of 'fn[[2]]'.")
            # cvfn <- list(sd(fn[[2]])/mfn[[1]], NULL)
            cvfn <- rep(sd(fn[[2]])/mean(fn[[2]]), 2)
          }
        }
      }

      if(!is.null(fitness)){
        if(sum(fitness < 0) != 0 | length(fitness) > 2) stop("if given, 'fitness' must be a vector with size <= 2 of positive numbers.")
        fitness0 <- NULL
      } else {
        if(missing(fitness0)) fitness0 <- 0
      }
      if(missing(mutations0)) mutations0 <- 0
      if(with.prob & missing(mutprob0)) mutprob0 <- 0

      if(sum(death < 0 | death >= 0.5) != 0 | length(death) > 2) stop("'death' must be a vector with size <= 2 of positive and < 0.5 numbers.")
      if(sum(plateff > 1) != 0 | length(plateff) > 2) stop("'plateff' must be a vector with size <= 2 of positive and <= 1 numbers.")
      if(sum(plateff < 1) != 0 & model == "H" & method == "ML"){
        warning("'plateff' can not be taking into account for ML method when 'model' is 'H': 'plateff' is set to 1.")
        plateff <- 1
      }
      if(sum(plateff < 1) != 0 & method == "P0"){
        warning("'plateff' can not be taking into account for P0 method: 'plateff' is set to 1.")
        plateff <- 1
      }
    }

    if(length(mutations0) > 1) {
      if(sum(mutations0 == mutations0[1]) == length(mutations0)) mutations0 <- mutations0[1]
      else stop("If you use a 'data.table', 'mutations0' can not have different values in a same class.")
    }
    if(mutations0 < 0) stop("'mutations0' must be a positive number.")

    if(!is.null(fitness0)){
      if(length(fitness0) > 1) {
        if(sum(fitness0 == fitness0[1]) == length(fitness0)) fitness0 <- fitness0[1]
        else stop("If you use a 'data.table', 'fitness0' can not have different values in a same class.")
      }
      if(fitness0 < 0) stop("'fitness0' must be a positive number.")
    }

    if(!is.null(mutprob0)){
      if(length(mutprob0) > 1) {
        if(sum(mutprob0 == mutprob0[1]) == length(mutprob0)) mutprob0 <- mutprob0[1]
        else stop("If you use a 'data.table', 'mutprob0' can not have different values in a same class.")
      }
      if(mutprob0 < 0 | mutprob0 >= 1) stop("if given, 'mutprob0' must be a positive and <= 1 number.")
      with.prob <- TRUE
      if(is.null(mfn))  mfn <- 1e9
      if(is.null(cvfn))  cvfn <- 0
    } else {
      if(with.prob){
        if(is.null(mfn))  mfn <- 1e9
        if(is.null(cvfn)) cvfn <- 0
        if(nsamples == 1) mutprob0 <- mutations0/mfn
        else mutprob0 <- 0
      }
    }

    if(length(conf.level) > 1) {
      if(sum(conf.level == conf.level[1]) == length(conf.level)) conf.level <- conf.level[1]
      else stop("If you use a 'data.table', 'conf.level' can not have different values in a same class.")
    }
    if(conf.level > 1 | conf.level < 0) stop("'conf.level' must be a positive and <= 1 numbers.")

    if(length(winsor) > 1) {
      if(sum(winsor == winsor[1]) == length(winsor)) winsor <- winsor[1]
      else stop("If you use a 'data.table', 'winsor' can not have different values in a same class.")
    }
    if(winsor < 0 | trunc(winsor) != winsor) stop("'winsor' must be a single positive integer.")

    if(model == "I"){
      if(!is.finite(muinf)){
        warning("'muinf' is infinite, 'model' is set to 'LD'.")
        model <- "LD"
      }
    } else {
      if(is.finite(muinf)) {
      	warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      	muinf <- +Inf
      }
    }


    H0 <- c(if(with.prob) mutprob0 else mutations0, fitness0)     # Vector of null hypothesises

    np <- length(H0)                         # Number of tested values

    if(missing(alternative)) alternative <- rep("two.sided", np)
    alternative <- match.arg(alternative, several.ok = TRUE)



    if(nsamples == 1){

      parameter <- NULL
      names <- character()

      if(np == 1){
        names(H0) <- if(with.prob) "mutation probability" else "mutation number"
        parameter <- c(fitness, death, plateff)
        names <- c("fitness", "death", "plateff")
        if(with.prob){
          if(is.list(mfn)) mfn <- unlist(mfn)
          if(is.list(cvfn)) cvfn <- unlist(cvfn)
          parameter <- c(parameter, mfn, cvfn)
          names <- c(names, "mfn", "cvfn")
        }
        if(is.finite(muinf)){
          parameter <- c(parameter, muinf)
          names <- c(names, "muinf")
        }
      }
        if(np == 2) {
        names(H0) <- c(if(with.prob) "mutation probability" else "mutation number", "fitness")
        parameter <- c(death, plateff)
        names <- c("death", "plateff")
        if(with.prob){
          if(is.list(mfn)) mfn <- unlist(mfn)
          if(is.list(cvfn)) cvfn <- unlist(cvfn)
          parameter <- c(parameter, mfn, cvfn)
          names <- c(names, "mfn", "cvfn")
        }
        if(model == "I"){
          parameter <- c(parameter, muinf)
          names <- c(names, "muinf")
        }
      }
      names(parameter) <- names
    }
    if(nsamples == 2){
      if(np == 1 & length(fitness) == 1) fitness <- c(fitness, fitness)
      if(length(death) == 1) death <- c(death, death)
      if(length(plateff) == 1) plateff <- c(plateff, plateff)
      if(length(muinf) == 1) muinf <- c(muinf, muinf)
      # if(length(mfn) == 1) mfn <- c(mfn, mfn)
      # if(length(cvfn) == 1) cvfn <- c(cvfn, cvfn)

      parameter <- NULL
      names <- character()


      if(np == 1){
        names(H0) <- if(with.prob) "mutprob difference" else "mutations difference"
        parameter <- cbind(fitness, death, plateff)
        names <- c("fitness", "death", "plateff")
        if(with.prob) {
          if(is.list(mfn)) mfn <- unlist(mfn)
          if(is.list(cvfn)) cvfn <- unlist(cvfn)
          parameter <- cbind(parameter, mfn, cvfn)
          names <- c(names, "mfn", "cvfn")
        }
        if(model == "I"){
          parameter <- cbind(parameter, muinf)
          names <- c(names, "muinf")
        }
      }
      if(np == 2){
        names(H0) <- c(if(with.prob) "mutprob difference" else "mutations difference", "fitness difference")
        parameter <- cbind(death, plateff)
        names <- c("death", "plateff")
        if(with.prob){
          if(is.list(mfn)) mfn <- unlist(mfn)
          if(is.list(cvfn)) cvfn <- unlist(cvfn)
          parameter <- cbind(parameter, mfn, cvfn)
          names <- c(names, "mfn", "cvfn")
        }
        if(model == "I"){
          parameter <- cbind(parameter, muinf)
          names <- c(names, "muinf")
        }
      }
      colnames(parameter) <- names
    }

    if(is.null(mfn))  mfn <- list(mfn)
    if(is.null(cvfn))  cvfn <- list(cvfn)
    if(is.null(fitness)) fitness <- list(fitness)

    # cat("length(mc) = ", length(mc), "\n")
    # cat("length(fn) = ", length(fn), "\n")
    # cat("mfn = ", mfn, "\n")
    # cat("mode(cvfn) = ", mode(cvfn), "\n")
    # cat("length(mfn) = ", length(mfn), "\n")
    # cat("length(cvfn) = ", length(cvfn), "\n")
    # cat("mfn = ", mfn[[1]], ", ", mfn[[2]], "\n")
    # cat("cvfn = ", cvfn[[1]], ", ", cvfn[[2]], "\n")

  # Estimates mean number of mutations alpha, given fitness parameter rho
    ests <- mapply(
            function(x, y, m, c, f, d, z, mui){
                    mutestim(mc = x, fn = y,
                            mfn = m, cvfn = c,
                            fitness = f, death = d, plateff = z,
                            model = model, muinf = mui,
                            method = method, winsor = winsor)
            }, mc, fn, mfn, cvfn, fitness, death, plateff, muinf)

    if(np == 1){   # If only mutations (or mutprob) is tested
      if(length(alternative) > 1) alternative <- alternative[1]
      if(nsamples == 2){    # If Two-sample test
        Tstat <- unlist(ests[1, ])    # Extract the estimate to compute statistic of the test
        sds <- unlist(ests[2, ])      # Standard deviation of the estimate
        ests <- Tstat                # Keep the estimate
        Tstat <- -diff(Tstat)      # Statistic of the test is build with difference of estimates
        sds <- sqrt(sum(sds^2))    # Standard deviation of the difference between estimates
      } else {
        Tstat <- unlist(ests[1, 1])    # Extract the estimate to compute statistic of the test
        sds <- unlist(ests[2, 1])      # Standard deviation of the estimate
        ests <- Tstat                # Keep the estimate
      }
    } else if(np == 2){    # If fitness is also tested
      if(length(alternative) == 1) alternative <- rep(alternative, 2)
      if(length(alternative) > 2) alternative <- alternative[c(1, 2)]
      if(nsamples == 2){               # If Two-sample test
        Tstat <- rbind(unlist(ests[1, ]), unlist(ests[3, ]))   # Extract estimates to compute statistic of the test
        sds <- rbind(unlist(ests[2, ]), unlist(ests[4, ]))     # Standard deviation of the estimates
        ests <- t(Tstat)                                      # Keep the estimates
        Tstat <- -apply(Tstat, 1, diff)                     # Statistic of the test is build with difference of estimates
        sds <- apply(sds, 1, function(s){sqrt(sum(s^2))})   # Standard deviation of the difference between estimates
      } else {
        Tstat <- c(unlist(ests[1, ]), unlist(ests[3, ]))   # Extract estimates to compute statistic of the test
        sds <- c(unlist(ests[2, ]), unlist(ests[4, ]))     # Standard deviation of the estimates
        ests <- Tstat                  # Keep the estimates
      }
    }

    cint <- mapply(
            function(e, s, alt){
              if(alt == "less"){                                        # Confidence interval(s) of the test
                c(-Inf, e+s*qnorm(conf.level))                           # with respect to the confidence level
              } else if(alt == "greater"){                              # and the alternative
                c(e-s*qnorm(conf.level), Inf)
              } else if(alt == "two.sided"){
                res <- e+s*c(-1, 1)*qnorm((1+conf.level)/2)
                if(nsamples == 1 & res[1] < 0) res[1] <- 0
                res
              }
            }, Tstat, sds, alternative)
    if(nsamples == 1) cint[1, cint[1, ] < 0] <-0

    Tstat <- (Tstat-H0)/sds                                      # Statistic(s) of the test

    pval <- mapply(
            function(alt, tstat){
              if(alt == "less"){                                        # p-value(s) of the test
                pnorm(tstat)                                            # with respect to the alternative
              } else if(alt == "greater"){
                pnorm(tstat, lower.tail = FALSE)
              } else if(alt == "two.sided"){
                2*pnorm(-abs(tstat))
              }
            }, alternative, Tstat)

    if(nsamples == 1){
      names(ests) <- c(if(with.prob) "mutation probability" else "mutation number", if(np == 2)"fitness")
    } else {
      if(np == 1) ests <- rbind(ests[1], ests[2])
      colnames(ests) <- c(if(with.prob) "mutation probability" else "mutation number", if(np == 2)"fitness")
    }
    # colnames(ests) <- names(H0)
    names(Tstat) <- names(H0)
    names(pval) <- names(H0)
    colnames(cint) <- names(H0)
    rownames(cint) <- c("bInf", "bSup")
    attr(cint, "conf.level") <- conf.level
  }
  # Build object with 'flantest' class
  flantest(Tstat = Tstat, estimate = ests, parameter = parameter, conf.int = cint, p.value = pval,
  null.value = H0, alternative = alternative, data.name = dname, model = model, method = method, nsamples = nsamples)
}



    ## --------------------------------------------------------------------------------------------------------------
    ############################ Density, distribution function, quantile function and random #######################
    ################################## generation for the mutants counts ############################################
    ## --------------------------------------------------------------------------------------------------------------

    # There are two cases:
    #  - If alpha is given, then returns a sample mcs of mutants counts with parameters:
    #   * mutations: mean number of mutation of mutation
    #   * fitness: fitness parameter
    #   * death: probability of death
    #   * dist: lifetimes distribution
    #
    #  - If alpha is not given, then compute first a sample fn of final counts with a Log-Normal distribution with mean mfn and
    #    coefficient of variation cvfn. Then compute a sample mcs of mutants counts, where the ith count has the following parameters:
    #   * mutprob*fni: mean number of mutation
    #                      (where fni is the corresponding of final count, and pi the probability of mutation)
    #   * fitness: fitness parameter
    #   * death: probability of death
    #   *dist: lifetimes distribution
    #
    #   Returns the two samples mcs and fn
    #
rflan <- function(n, mutations = 1., mutprob = NULL, fitness = 1., death = 0., plateff = 1.,
        dist = list(name = "lnorm", meanlog = -0.3795851, sdlog = 0.3016223),
        distfn = NULL,
        mfn = 1.e9, cvfn = 0,
        muih = list(mu = NULL, muinv0 = NULL),
        ...) {

    if(n <= 0) stop("'n' must be a positive integer.")

    if(length(n) > 1) n <- length(n)

    if(mutations < 0 | length(mutations) > 1) stop("'mutations' must be a single positive number")

    if(fitness < 0 | length(fitness) > 1) stop("'fitness' must be a single positive number")

    if(mfn < 0 | length(mfn) > 1) stop("'mfn' must be a single positive number")

    if(cvfn < 0 | length(cvfn) > 1) stop("'cvfn' must be a single positive number")

    if((!is.list(muih)) | (length(muih) > 2)) stop("'muih' must be a list of a function which can be followed by its inverse.")

    if(!is.list(dist)){
      if(dist == "exp") {
        dist <- list("exp")
      } else if(dist == "dirac") {
        dist <- list("dirac")
      } else stop("'dist' must be a list of a character chain followed by its arguments if required")
    }


    names(dist)[1] <- "name"
    # if(dist$name == "exp"){
    #   names(dist)[2] <- "rate"
    # } else if(dist$name == "dirac"){
    #   names(dist)[2] <- "location"
    # } else
    if(dist$name == "lnorm"){
      names(dist)[2] <- "meanlog"
      names(dist)[3] <- "sdlog"
    } else if(dist$name == "gamma"){
      names(dist)[2] <- "shape"
      names(dist)[3] <- "scale"
    } else if(dist$name != "exp" & dist$name != "dirac") stop("'dist[[1]]' must be a character chain among 'exp', 'dirac', 'lnorm', 'gamma'.")

    if(cvfn > 0){
      if(!is.null(distfn)){
       if(distfn != "lnorm" & distfn != "gamma") stop("'distfn' must be a character chain among 'lnorm', 'gamma'.")
     } else distfn <- "lnorm"
    } else {
      if(!is.null(distfn)) {
        warning("'cvfn' is zero: 'distfn' is ignored.")
      }
      distfn <- "null"
    }

    if(is.null(muih[[1]])){
      muih <- NULL
    } else {
      match.fun(muih[[1]])
      mu <- muih[[1]]
      MU <- function(u, v) mu(u, v, ...)

      muinf <- MU(0, +Inf)
      if(is.na(muinf)) stop("limit of 'muih[[1]]' is not defined (NaN)")
      if(!is.finite(muinf)){
      	dist <- list(name = "exp", 1)
      	muih <- NULL
      }
      if((length(muih) == 1)){
        muinv0 <- function(u) {
          x <- seq(MU(0, u), muinf, by = MU(0, u)/muinf)
          fx <- sapply(x, MU)
          ind <- which(fx >= u)
          x[min(ind)]
        }
        muih <- list(mu = MU, muinv0 = muinv0)
      } else if(length(muih) == 2){
        match.fun(muih[[2]])
        muinv0 <- muih[[2]]
        MUINV0 <- function(u) muinv0(u, ...)
        muih <- list(mu = MU, muinv0 = MUINV0)
      }
    }


    if(death < 0 | death >= 0.5 | length(death) > 1) stop("'death' must be a single positive and < 0.5 number ")

    if(plateff < 0 | plateff > 1 | length(plateff) > 1) stop("'plateff' must be a single positive and <= 1 number.")

    if(!is.null(mutprob)){
      if(mutprob < 0 | mutprob > 1 | length(mutprob) > 1) stop("'mutprob' must be a single positive and <= 1 number")
      mutations <- mutprob*(if(cvfn == 0) mfn else 1)              # Sample with mutprob instead of mutations if cvfn > 0
    } else {
      if(cvfn > 0) mutations <- mutations/mfn    # Default value of mutprob if missing and cvfn > 0
    }

    if(!is.null(muih)){
      flan.sim <- new(FlanIhSim, list(
      	mutations = mutations,
      	fitness = fitness,
      	death = death,
      	muih = muih,
        distfn = distfn,
      	mfn = mfn,
      	cvfn = cvfn
      ))
    } else {
      if(dist$name == "lnorm" | dist$name == "gamma")  dist <- adjust.rate(dist, 1, death)
      flan.sim <- new(FlanSim, list(
        mutations = mutations,
        fitness = fitness,
        death = death,
        dist = dist,
        distfn = distfn,
        mfn = mfn,
        cvfn = cvfn
      ))
    }
    output <- flan.sim$rflan(n)


    mc <- output$mc
    indNA <- which(mc < 0)
    if(length(indNA) > 0){
      warning("Production of NaNs in geometric sampling.")
      output$mc[indNA] <- NaN
    }

    if(cvfn == 0) output$fn <- rep(mfn, n)
    if(plateff < 1) output$mc <- sapply(output$mc, function(mc) if(mc > 0) rbinom(1, size = mc, prob = plateff) else 0)

    output
}

  # Quantile function of the mutants count with parameter
  #   mutations: mean number of mutations
  #   fitness: fitness parameter
  #   death: death probability
  #   model: lifetimes distribution model: exponentialy distributed lifetimes ("LD") or constant lifetimes ("H")

qflan <- function(p, mutations = 1., fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, lower.tail = TRUE){

  if((sum(p < 0)+sum(p > 1)) != 0){
    stop("'p' must be a vector of positive and <= 1 numbers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  if(missing(model)){model = "LD"}
  model <- match.arg(model)

  if(plateff < 0 | plateff > 1| length(plateff) > 1){
    stop("'plateff' must be a single positive and <= 1 number.")
  }

  if(model == "H" & plateff < 1){
    # stop("If 'model' is 'H', 'plateff' can not be < 1.")
    warning("Probabilities are not available when 'model' is 'H' and 'plateff' < 1. 'plateff' is set to 1")
    plateff <- 1
  }

  m <- 100
  P <- pflan(m = 0:m, mutations = mutations, fitness = fitness, death = death, plateff = plateff, model = model, lower.tail = lower.tail)
  sapply(p, function(pp){
    if(pp == 1) k <- Inf
    if (lower.tail & pp <= P[1]) k <- 0
    if (!lower.tail & pp >= P[1]) k <- 0
    else {
      k <- if(lower.tail) max(which(P<pp)) else max(which(P>pp))          # quantile if in the actual table
      while (k>= m){                   # if p not yet in the table
	m2 <- 2*m                       # double the table
	P2 <- pflan(m = m:m2, mutations = mutations, fitness = fitness, death = death, plateff = plateff, model = model, lower.tail = lower.tail)      # truncated distribution function
	if (lower.tail & pp <= P2[1]) k
	if (!lower.tail & pp >= P2[1]) k
	else {
	  k <- k-1+if(lower.tail) max(which(P2<pp)) else max(which(P2>pp))           # quantile if in the table
	  m <- m2
	}
      }                              # end while
    }                                  # end if
    k
  })

}

  # Distribution function of the mutants count with parameter
  #   mutations: mean number of mutations
  #   fitness: fitness parameter
  #   death: death probability
  #   model: lifetimes distribution model: exponentialy distributed lifetimes ("LD") or constant lifetimes ("H")


pflan <- function(m, mutations = 1., fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, lower.tail = TRUE){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  if(missing(model)){model = "LD"}
  model <- match.arg(model)

  if(plateff < 0 | plateff > 1| length(plateff) > 1){
    stop("'plateff' must be a single positive and <= 1 number.")
  }
  if(model == "I"){
    if(!is.finite(muinf)){
      warning("'muinf' is infinite, 'model' is set to 'LD'.")
      model <- "LD"
    }
    # if(plateff < 1) model <- "Ipef"
  } else {
    if(is.finite(muinf)) {
      warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      muinf <- Inf
    }
  }

  if(model == "H" & plateff < 1){
    # stop("If 'model' is 'H', 'plateff' can not be < 1.")
    warning("Probabilities are not available when 'model' is 'H' and 'plateff' < 1. 'plateff' is set to 1")
    plateff <- 1
  }
  # if(model == "LD" & plateff < 1) model <- "LDpef"


  flan.mutmodel <- new(FlanMutMod, list(
    mutations = mutations,
    fitness = fitness,
    death = death,
    plateff = plateff,
    # integrands = integrands,
    model = model,
    muinf = muinf,
    mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
  ))

  M <- max(m)
  output <- flan.mutmodel$pflan(M, lower.tail)

  output[m+1]
}


  # Density function of the mutants count with parameter
  #   mutations: mean number of mutations
  #   fitness: fitness parameter
  #   death: death probability
  #   model: lifetimes distribution model: exponentialy distributed lifetimes ("LD") or constant lifetimes ("H")


dflan <- function(m, mutations = 1., fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }

  if(missing(model)){model = "LD"}
  model <- match.arg(model)

  if(plateff < 0 | plateff > 1| length(plateff) > 1){
    stop("'plateff' must be a single positive and <= 1 number.")
  }

  if(model == "I"){
    if(!is.finite(muinf)){
      warning("'muinf' is infinite, 'model' is set to 'LD'.")
      model <- "LD"
    }
    # if(plateff < 1) model <- "Ipef"
  } else {
    if(is.finite(muinf)) {
      warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      muinf <- Inf
    }
  }
  if(model == "H" & plateff < 1){
    warning("Probabilities are not available when 'model' is 'H' and 'plateff' < 1. 'plateff' is set to 1")
    plateff <- 1
  }
  # if(model == "LD" & plateff < 1) model <- "LDpef"


  flan.mutmodel <- new(FlanMutMod, list(
    mutations = mutations,
    fitness = fitness,
    death = death,
    plateff = plateff,
    # integrands = integrands,
    model = model,
    muinf = muinf,
    mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
  ))

  M <- max(m)
  output <- flan.mutmodel$dflan(M)

  output[m+1]
}



    ## --------------------------------------------------------------------------------------------------------------
    ############################################## Hidden functions #################################################
    ## --------------------------------------------------------------------------------------------------------------


## Sampling
adjust.rate <- function(dist, fitness = 1., death = 0.){
#   rescales the parameter(s) of distribution dist to obtain growh rate equal to fitness.
#   The distribution is given and returned as a list of a character chain
#   followed by the list of parameters.

  adist <- dist				# adjusted distribution
  m <- 2*(1-death)			# mean offspring number
  switch(dist[[1]],
    dirac = {                              # Dirac distribution
      adist <- list(name = "dirac", location = log(m)/fitness) # new location parameter
	  },                              # end Dirac distribution
    exp = {                                # Exponential distribution
    	adist <- list(name = "exp", rate = fitness/(m-1))
	  },
    gamma = {				# Gamma distribution
      a <- dist$shape			# according to the shape parameter
      adist$scale <- (m^(1/a)-1)/fitness	# adjust scale parameter
    },
    lnorm = {					# Log-normal distribution
      a <- dist$sdlog
      l <- dist$meanlog
      lap <- function(s){               # Laplace transform
	f <- function(x){exp(-(s*exp(a*x*sqrt(2)+l)+x^2))}

	I <- {integrate(f,
	      lower = -Inf, upper = Inf, rel.tol = .Machine$double.eps^0.5, subdivisions = 1.e3L)}

	I <- I$value                      # get the value
  	I <- I/sqrt(pi)			# rescale
        I-1/m                     # return shifted Laplace transform
      }
      lwb <- 1                          # lower bound of search interval
      while(lap(lwb)<0){lwb<-lwb/2}     # adjust lower bound
      upb <- 1                          # upperbound of search interval
      while(lap(upb)>0){upb<-upb+1}     # adjust upperbound
      s <- uniroot(lap, c(lwb, upb))$root       # find root of lap
      adist$meanlog <- log(s/fitness)+adist$meanlog	# adjust scale parameter

    }
  )
  adist
}

## Estimation

	      #//////////////////////////////////ML method//////////////////////////////////#

# Returns the ML estimate of mean number of mutation for a sample mc, given the fitness and death
# If mfn or cvfn are non-empty, returns the estimate of mutation probability instead of the mean number, after decreasing the bias
# induced if cvfn > 0
MutationMLOptimization <- function(mc, mfn = NULL, cvfn = NULL, fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, winsor = 2000){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)


  # Initialization
  est <- MutationGFEstimation(mc, fitness = fitness, death = death, plateff = plateff, model = model, muinf = muinf, init = TRUE)

  a.est <- est$mutations

  # Winsorization
  if(max(mc) > winsor) mc[mc > winsor] = winsor

  if(dflan(m = max(mc), mutations = a.est, fitness = fitness, death = death, plateff = plateff, model = model, muinf = muinf) == 0) {
    warning("Infinite log-likelihood: returns GF estimates (try a larger value for 'winsor')." )
    return(list(mutations = a.est, sd.mutations = est$sd.mutations))
  }

  dC <- dclone(m = 0:max(mc), fitness = fitness, death = death, plateff = plateff, model = model)

  ll <- function(a){
    p <- log(deduce.dflan(m = mc, mutations = a, fitness = fitness, death = death, clone = dC))
    -sum(p)
  }
  # gradient of log-likelihood of the sample
  dll <- function(a){
    p <- deduce.dflanda(m = mc, mutations = a, fitness = fitness, death = death, clone = dC)
    res <- p$dQ_da/p$Q
    -sum(res)
  }


  lower = 0.1*a.est
  upper = 10*a.est

  # a.est <- lbfgsb3c(par = a.est, fn = ll, gr = dll, lower = lower, upper = upper, control = list(trace = -1))$par
  a.est <- optim(par = a.est, fn = ll, gr = dll, method="L-BFGS-B", lower = lower, upper = upper, control = list(trace = 0))$par

  dldd <- deduce.dflanda(m = mc, mutations = a.est, fitness = fitness, death = death, clone = dC)

  I <- sum((dldd$dQ_da)^2/(dldd$Q)^2)

  sda <- sqrt(1/I)

  if(!is.null(mfn)) {
    if(cvfn > 0){
#       z4 = .tunings$z4
      z4 <- 0.55
  #     if(model == "LD") {
	# integrands <- list(CLONE_PGF = function(x, rho, delta) {x^rho/(1+x*delta)})
  #     } else integrands <- NULL
      Mutmodel <- new(FlanMutMod, list(
		mutations = a.est,
		fitness = fitness,
		death = death,
		# integrands = integrands,
		model = model,
    muinf = muinf,
    plateff = NULL, mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
	      ))
      pm.est <- Mutmodel$unbias.mutprob(sda, z4, mfn, cvfn)
    } else pm.est <- list(mutprob = a.est/mfn, sd.mutprob = sda/mfn)
    pm.est
  } else list(mutations = a.est, sd.mutations = sda)


}

# Returns the ML estimate of mutation probability for a sample of couple (mc, fn), given the fitness and death
MutationProbabilityMLOptimization <- function(mc, fn, fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, winsor = 2000){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  mfn <- mean(fn)
  cvfn <- sd(fn)/mfn

  # Initialization
  est <- MutationGFEstimation(mc = mc, mfn = mfn, cvfn = cvfn, fitness = fitness, death = death, plateff = plateff, model = model, muinf = muinf, init = TRUE)

  pm.est <- est$mutprob*mfn
  # Winsorization
  if(max(mc) > winsor) mc[mc > winsor] = winsor

  if(dflan(m = max(mc), mutations = pm.est, fitness = fitness, death = death, plateff = plateff, model = model, muinf = muinf) == 0) {
    warning("Infinite log-likelihood: returns GF estimates (try a larger value for 'winsor')." )
    return(list(mutprob = pm.est/mfn, sd.mutprob = est$sd.mutprob))
  }

  dC <- dclone(m = 0:max(mc), fitness = fitness, death = death, plateff = plateff, model = model, muinf = muinf)

  ll <- function(pm){
    p <- mapply(function(x, y){
      y <- y/mfn
      log(deduce.dflan(m = x, mutations = pm*y, fitness = fitness, death = death, clone = dC))
    }, mc, fn)
    -sum(p)
  }

  # gradient of log-likelihood of the sample
  dll <- function(pm){
    res <- mapply(function(x, y){
      y <- y/mfn
      p<-deduce.dflanda(m = x, mutations = pm*y, fitness = fitness, death = death, clone = dC)
      p$dQ_da*y/p$Q
    }, mc, fn)
    -sum(res)
  }

  lower = 0.1*pm.est
  upper = 10*pm.est

  # pm.est <- lbfgsb3(par = pm.est, fn = ll, gr = dll, lower = lower, upper = upper, control = list(trace = -1))$par/mfn
  pm.est <- optim(par = pm.est, fn = ll, gr = dll, method="L-BFGS-B", lower = lower, upper = upper, control = list(trace = 0))$par/mfn

  dldd <- mapply(function(x, y){
	    deduce.dflanda(m = x, mutations = pm.est*y, fitness = fitness, death = death, clone = dC)
	  }, mc, fn)
  dldd <- list(Q = unlist(dldd[1, ]), dQ_da = unlist(dldd[2, ]))
  I <- sum((dldd$dQ_da*fn)^2/dldd$Q^2)            # Fisher information

  sdpm <- sqrt(1/I)

  list(mutprob = pm.est, sd.mutprob = sdpm)

}

# Returns the ML estimates of mean number of mutation and fitness for a sample mc, given the death
# If mfn or cvfn are non-empty, returns the estimate of mutation probability instead of the mean number, after decreasing the bias
# induced if cvfn > 0
MutationFitnessMLOptimization <- function(mc, mfn = NULL, cvfn = NULL, death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, winsor = 2000){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  # Initialization
  est <- MutationFitnessGFEstimation(mc = mc, death = death, plateff = plateff, model = model, muinf = muinf, init = TRUE)
  if(!est$succeeds) warning("Initialization of 'fitness' with 'GF'-method has failed.")
  a.est <- est$mutations ; r.est <- est$fitness

  # Winsorization
  if(max(mc) > winsor) mc[mc > winsor] = winsor

  if(dflan(m = max(mc), mutations = a.est, fitness = r.est, death = death, plateff = plateff, model = model, muinf = muinf) == 0) {
    warning("Infinite log-likelihood: returns GF estimates (try a larger value for 'winsor')." )
    return(list(mutations = a.est, sd.mutations = est$sd.mutations, fitness = r.est, sd.fitness = est$sd.fitness))
  }

  ll <- function(param){
    a <- param[1]
    r <- param[2]
    p <- dflan(m = mc, mutations = a, fitness = r, death = death, plateff = plateff, model = model, muinf = muinf)
    p <- log(p)
    -sum(p)
  }

  # gradient of log-likelihood of the sample
  dll <- function(param){
    a = param[1]
    r = param[2]
    p <- dflan.grad(m = mc, mutations = a, fitness = r, death = death, plateff = plateff, model = model, muinf = muinf, dalpha = TRUE, drho = TRUE)
    res <- rbind(p$dQ_da/p$Q, p$dQ_dr/p$Q)
    -apply(res, 1, sum)
  }

  lower = 0.1*c(a.est, r.est)
  upper = 10*c(a.est, r.est)

#   if(!est$succeeds) {
#     lower[2] <- 10*r.est ;
#     upper[2] <- 500*r.est
#   }
#
  # est <- lbfgsb3(par = c(a.est, r.est), fn = ll, lower = lower, upper = upper, control = list(trace = -1))$par
  est <- optim(par = c(a.est, r.est), fn = ll, gr = dll, method="L-BFGS-B", lower = lower, upper = upper, control = list(trace = 0))$par

  a.est <- est[1]					# Update alpha estimate
  r.est <- est[2]					# Update rho estimate

  dldd <- dflan.grad(m = mc, mutations = a.est, fitness = r.est, death = death, plateff = plateff, model = model, muinf = muinf, dalpha = TRUE, drho = TRUE)

  p2 <- (dldd$Q)^2
  dpa <- dldd$dQ_da
  dpr <- dldd$dQ_dr

  Ia <- sum((dpa^2/p2))
  Ir <- sum(dpr^2/p2)
  Iar <- sum(dpa*dpr/p2)
# #
  det <- Ia*Ir-Iar^2
  sda <- sqrt(Ir/det)
  sdr <- sqrt(Ia/det)

  if(!is.null(mfn)) {
    if(cvfn > 0){
      z4 <- 0.55
      # if(model == "LD") {
	# integrands <- list(CLONE_PGF = function(x, rho, delta) {x^rho/(1+x*delta)})
  #     } else integrands <- NULL
      Mutmodel <- new(FlanMutMod, list(
		mutations = a.est,
		fitness = r.est,
		death = death,
		# integrands = integrands,
		model = model,
    muinf = muinf,
    plateff = NULL, mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
	      ))
      pm.est <- Mutmodel$unbias.mutprob(sda, z4, mfn, cvfn)
    } else pm.est <- list(mutprob = a.est/mfn, sd.mutprob = sda/mfn)

    c(pm.est, fitness = r.est, sd.fitness = sdr)
  } else list(mutations = a.est, sd.mutations = sda, fitness = r.est, sd.fitness = sdr)


}


# Returns the ML estimates of mutation probability and fitness for a sample of couple (mc, fn), given the death
MutationProbabilityFitnessMLOptimization <- function(mc, fn, death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, winsor = 2000){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  mfn <- mean(fn)
  cvfn <- sd(fn)/mfn


  # Initialization
  est <- MutationFitnessGFEstimation(mc, mfn = mfn, cvfn = cvfn, death = death, plateff = plateff, model = model, muinf = muinf, init = TRUE)
  if(!est$succeeds) warning("Initialization of 'fitness' with 'GF'-method has failed.")
  pm.est <- est$mutprob*mfn ; r.est <- est$fitness

  # Winsorization
  if(max(mc) > winsor) mc[mc > winsor] = winsor

  if(dflan(m = max(mc), mutations = pm.est, fitness = r.est, death = death, plateff = plateff, model = model, muinf = muinf) == 0) {
    warning("Infinite log-likelihood: returns GF estimates (try a larger value for 'winsor')." )
    return(list(mutprob = pm.est/mfn, sd.mutations = est$sd.mutprob, fitness = r.est, sd.fitness = est$sd.fitness))
  }

  ll <- function(param){
    pm <- param[1]
    r <- param[2]
    p <- mapply(function(x, y){
      y <- y/mfn
      log(dflan(m = x, mutations = pm*y, fitness = r, death = death, plateff = plateff, model = model, muinf = muinf))
    }, mc, fn)
    -sum(p)
  }

  # gradient of log-likelihood of the sample
  dll <- function(param){
    pm = param[1]
    r = param[2]
    res <- mapply(function(x, y){
      y <- y/mfn
      p<-dflan.grad(m = x, mutations = pm*y, fitness = r, death = death, plateff = plateff, model = model, muinf = muinf, dalpha = TRUE, drho = TRUE)
      rbind(p$dQ_da*y/p$Q, p$dQ_dr/p$Q)
    }, mc, fn)
    -apply(res, 1, sum)
  }


  lower = 0.1*c(pm.est, r.est)
  upper = 10*c(pm.est, r.est)

#   if(!est$succeeds) {
#     lower[2] <- 10*r.est ;
#     upper[2] <- 500*r.est
#   }

  # est <- lbfgsb3(par = c(pm.est, r.est), fn = ll, gr = dll, lower = lower, upper = upper, control = list(trace = -1))$par
  est <- optim(par = c(pm.est, r.est), fn = ll, gr = dll, method="L-BFGS-B", lower = lower, upper = upper, control = list(trace = 0))$par
  pm.est = est[1]/mfn					# Update alpha estimate
  r.est = est[2]					# Update rho estimate

  dldd <- mapply(function(x, y){
	    dflan.grad(m = x, mutations = pm.est*y, fitness = r.est, death = death, plateff = plateff, model = model, muinf = muinf, dalpha = TRUE, drho = TRUE)
	  }, mc, fn)

  p2 <- (unlist(dldd[1, ]))^2
  dppm <- unlist(dldd[2, ])*fn
  dpr <- unlist(dldd[3, ])

  Ipm <- sum((dppm^2/p2))
  Ir <- sum(dpr^2/p2)
  Ipmr <- sum(dppm*dpr/p2)
# #
  det <- Ipm*Ir-Ipmr^2
  sdpm <- sqrt(Ir/det)
  sdr <- sqrt(Ipm/det)

  list(mutprob = pm.est, sd.mutprob = sdpm, fitness = r.est, sd.fitness = sdr)

}


	      #//////////////////////////////////P0 method//////////////////////////////////#

# Returns the P0 estimate of mean number of mutations for a sample of couple mc, given the death
MutationsNumberP0Estimation <- function(mc, death = 0., plateff = 1.){
# TODO: simulation study for P0 with plateff
# Return the estimate of alpha for a sample mc
# by p0-method under cell deaths with probability delta
# # when delta is known.
  if(death > 0 | plateff < 1){
    dstar <- death/(1-death)			# extinction probability
    if(plateff < 1-dstar) {
      warning("'plateff < 1-death/(1-death)': Fixed point method cannot be applied. 'plateff' is set to 1.")
      plateff <- 1
    }
    epgf <- function(z) mean(z^mc)		# empirical PGF

    a <- -log(epgf((dstar-(1-plateff))/plateff))/(1-dstar)	#  estimate alpha

    sda <- (1-dstar)^(-2)*(epgf(((dstar-(1-plateff))/plateff)^2)/(epgf(dstar-(1-plateff)/plateff)^2)-1)	# variance of alpha's estimate

    sda <- sqrt(sda/length(mc))
  } else {
    p0 <- mean(mc == 0)
    a <- -log(p0)
    sda <- sqrt((1/p0-1)/length(mc))
  }

  list(mutations = a, sd.mutations = sda)

}


# Returns the P0 estimate of mean number of mutations for a sample of couple mc, given the death
# If mfn or cvfn are non-empty, returns the estimate of the mutation probability instaed of the mean number
# after decreasing the induced bias if cvfn > 0
MutationP0Estimation <- function(mc, fn = NULL, mfn = NULL, cvfn = NULL, death = 0.){

# Return the estimate of alpha for a sample mc
# by p0-method under cell deaths with probability delta
# # when delta is known.

  if(is.null(fn) | death > 0) {
    a <- MutationsNumberP0Estimation(mc = mc, death = death)

    sda <- a$sd.mutations
    a <- a$mutations

    if(!is.null(fn)){
      mfn <- mean(fn)
      cvfn <- sd(fn)/mfn
    }
    if(!is.null(mfn)){
      pm <- a/mfn
      sdpm <- sda/mfn
      if(cvfn != 0){
        umds <- 1-death/(1-death)			# extinction probability
        temp <- a*(cvfn^2)*umds         # relative bias on mutprob
        pm <- pm*(1+temp/2)		# unbiased estimator of mutrate
        sdpm <- sdpm*(1+temp)		# standard deviation of unbiased estimator
      }
      list(mutprob = pm, sd.mutprob = sdpm)
    } else list(mutations = a, sd.mutations = sda)
  } else MutationProbabilityP0MLOptimization(mc, fn)
}

# Returns the P0 estimate of mean number of mutations and fitness for a sample of couple mc, given the death
# If mfn or cvfn are non-empty, returns the estimate of the mutation probability instaed of the mean number
# after decreasing the induced bias if cvfn > 0
MutationFitnessP0Estimation <- function(mc, fn = NULL, mfn = NULL, cvfn = NULL, death = 0., model = c("LD", "H", "I"), muinf = +Inf, winsor = 2000){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  if(!is.null(fn) & death == 0){

    pm <- MutationProbabilityP0MLOptimization(mc, fn)
    sdpm <- pm$sd.mutprob
    pm <- pm$mutprob
    output <- list(mutprob = pm, sd.mutprob = sdpm)
    r <- FitnessP0Optimization(mc = mc, fn = fn, mut = pm, death = death, model = model, muinf = muinf, winsor = winsor)
  } else {
    a <- MutationsNumberP0Estimation(mc, death)

    sda <- a$sd.mutations
    a <- a$mutations

    if(!is.null(mfn)){
      pm <- a/mfn
      sdpm <- sda/mfn
      if(cvfn != 0){
        umds <- 1-death/(1-death)			# extinction probability
        temp <- a*(cvfn^2)*umds         # relative bias on mutprob
        pm <- pm*(1+temp/2)		# unbiased estimator of mutrate
        sdpm <- sdpm*(1+temp)		# standard deviation of unbiased estimator
      }
      output <- list(mutprob = pm, sd.mutprob = sdpm)
    } else output <- list(mutations = a, sd.mutations = sda)
    r <- FitnessP0Optimization(mc = mc, mut = a, death = death, model = model, winsor = winsor, muinf = muinf)
  }
  c(output, fitness = r$fitness, sd.fitness = r$sd.fitness)

}

# Returns the P0ML estimate of mutation probability for a sample of couple (mc, fn), given the death
MutationProbabilityP0MLOptimization <- function(mc, fn){

  mfn <- mean(fn)
  cvfn <- sd(fn)/mfn

  # Initialization
  est <- MutationP0Estimation(mc = mc, mfn = mfn, cvfn = cvfn, death = 0)

  pm.est <- est$mutprob*mfn

  X <- as.numeric(mc == 0)

  ll <- function(pm){
    p <- mapply(function(x, y){
      y <- y/mfn
      -pm*x*y+(1-x)*log(1-exp(-pm*y))
    }, X, fn)
    -sum(p)
  }

  # gradient of log-likelihood of the sample
  dll <- function(pm){
    res <- mapply(function(x, y){
      y <- y/mfn
      -x*y+(1-x)*y/(exp(pm*y)-1)
    }, mc, fn)
    -sum(res)
  }

  lower = 0.1*pm.est
  upper = 10*pm.est

  # pm.est <- lbfgsb3(par = pm.est, fn = ll, gr = dll, lower = lower, upper = upper, control = list(trace = -1))$par/mfn
  pm.est <- optim(par = pm.est, fn = ll, gr = dll, method="L-BFGS-B", lower = lower, upper = upper, control = list(trace = 0))$par/mfn

  dldd <- mapply(function(x, y){
	   -x*y+(1-x)*y/(exp(pm.est*y)-1)
  }, mc, fn)
  I <- sum(dldd^2)
  sdpm <- sqrt(1/I)

  list(mutprob = pm.est, sd.mutprob = sdpm)

}



# Returns the ML estimate of fitness for a sample of couple mc, given the mean number of mutations (or mutation probability) and death
# If fn is non-empty,
FitnessP0Optimization <- function(mc, fn = NULL, mut, death = 0., model = c("LD", "H", "I"), muinf = +Inf, winsor = 2000){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  # Initialization
  est <- FitnessGFEstimation(mc, death = death, plateff = 1., model = model, muinf = muinf)
  if(!est$succeeds) warning("Initialization of 'fitness' with 'GF'-method has failed.")
  r.est <- est$fitness

  # Winsorization
  if(max(mc) > winsor){ mc[which(mc>winsor)] = winsor}

  if(is.null(fn)){
    if(missing(mut)) mut <- 1

  if(dflan(m = max(mc), mutations = mut, fitness = r.est, death = death, model = model, muinf = muinf) == 0) {
    warning("Infinite log-likelihood: returns GF estimates (try a larger value for 'winsor')." )
    return(list(fitness = r.est, sd.fitness = est$sd.fitness))
  }

    ll <- function(r){
      p <- log(dflan(m = mc, mutations = mut, fitness = r, death = death, model = model, muinf = muinf))
      -sum(p)
    }

    # gradient of log-likelihood of the sample
    dll <- function(r){
      p <- dflan.grad(m = mc, mutations = mut, fitness = r, death = death, model = model, muinf = muinf, dalpha = FALSE, drho = TRUE)
      res <- p$dQ_dr/p$Q
      -sum(res)
    }
  } else {
    if(missing(mut)) mut <- 1e-9

  if(dflan(m = max(mc), mutations = mut*mean(fn), fitness = r.est, death = death, model = model) == 0) {
    warning("Infinite log-likelihood: returns GF estimates (try a larger value for 'winsor')." )
    return(list(fitness = r.est, sd.fitness = est$sd.fitness))
  }

    ll <- function(r){
      p <- mapply(function(x, y){
	    log(dflan(m = x, mutations = mut*y, fitness = r, death = death, model = model, muinf = muinf))
	  }, mc, fn)
      -sum(p)
    }

    # gradient of log-likelihood of the sample
    dll <- function(r){
      res <- mapply(function(x, y){
	      p<-dflan.grad(m = x, mutations = mut*y, fitness = r, death = death, model = model, muinf = muinf, dalpha = FALSE, drho = TRUE)
	      p$dQ_dr/p$Q
	    }, mc, fn)
      -sum(res)
    }

  }


  lower = 0.1*r.est
  upper = 10*r.est

  # r.est <- lbfgsb3(par = r.est, fn = ll, gr = dll, lower = lower, upper = upper, control = list(trace = -1))$par
  r.est <- optim(par = r.est, fn = ll, gr = dll, method="L-BFGS-B", lower = lower, upper = upper, control = list(trace = 0))$par
#   est2 <- optim(par = r.est, fn = ll, gr = dll, method = 'L-BFGS-B', lower = lower, upper = upper, hessian = TRUE)

  if(is.null(fn)) dldd <- dflan.grad(m = mc, mutations = mut, fitness = r.est, death = death, model = model, muinf = muinf, dalpha = FALSE, drho = TRUE)

  else {
    dldd <- mapply(function(x, y){
		  dflan.grad(m = x, mutations = mut*y, fitness = r.est, death = death, model = model, muinf = muinf, dalpha = FALSE, drho = TRUE)
	       }, mc, fn)
    dldd <- list(Q = unlist(dldd[1, ]), dQ_dr = unlist(dldd[2, ]))
  }

  I <- sum((dldd$dQ_dr)^2/(dldd$Q)^2)

  sdr <- sqrt(1/I)
  list(fitness = r.est, sd.fitness = sdr)
}

	      #//////////////////////////////////GF method//////////////////////////////////#

MutationGFEstimation <- function(mc, mfn = NULL, cvfn = NULL, fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, init = FALSE){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  q <- 0.1
  b <- quantile(mc, q, names = FALSE)+1           # scaling factor

  # if(model == "LD") {
  #   integrands <- list(CLONE_PGF = function(x, rho, delta) {x^rho/(1+x*delta)})
  # } else integrands <- NULL
  Mutmodel <- new(FlanMutMod, list(
    mc = mc,
    mfn = mfn, cvfn = cvfn,
    fitness = fitness, death = death, plateff = plateff,
    # integrands = integrands,
    model = model,
    muinf = muinf,
    scale = b,
    mutations = NULL
  ))

  Mutmodel$MutationGFEstimation(init)

}

FitnessGFEstimation <- function(mc, death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)

  q <- 0.1
  b <- quantile(mc, q, names = FALSE)+1           # scaling factor

  z1 <- 0.1                            # lower value
  z2 <- 0.9                            # higher value

  z1<-z1^(1/b); z2<-z2^(1/b) ;            # rescale variables


  g1 <- mean(z1^mc)                   # empirical generating function at z1
  g2 <- mean(z2^mc)                   # empirical generating function at z2
  y <- log(g1)/log(g2)                   # get ratio of logs

  # if(model == "LD") clone = new(FlanExpClone, list(death = death, integrands = list(CLONE_PGF = function(x, rho, delta) {x^rho/(1+x*delta)})))
  if(model == "LD") clone = new(FlanExpClone, list(death = death, fitness = NULL, plateff = NULL))
  else if(model == "H") clone = new(FlanDirClone, list(death = death, fitness = NULL))
  else clone=new(FlanInhClone, list(fitness = NULL, death = death, muinf = muinf, plateff = NULL))

  if(plateff < 1) {
    ump <- 1-plateff
    z1 <- ump+plateff*z1
    z2 <- ump+plateff*z2
  }

  f <- function(r){
    PGF <- clone$pgf2(r, c(z1, z2))
    ((1-PGF[1])/(1-PGF[2]))-y
  }

  binf <- 0.01
  bsup <- 100

  GFrho.succeeds <- TRUE
  if(f(binf)*f(bsup)>0) {
    rho <- 1
    GFrho.succeeds <- FALSE
  } else rho <- uniroot(f, interval = c(binf, bsup), tol = 1.e-6, maxiter = 1.e3)$root

  list(fitness = rho, succeeds = GFrho.succeeds)

}


MutationFitnessGFEstimation <- function(mc, mfn = NULL, cvfn = NULL, death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, init = FALSE){

  if(missing(model)) model <- "LD"
  model <- match.arg(model)


  q <- 0.1
  b <- quantile(mc, q, names = FALSE)+1           # scaling factor

  z1 <- 0.1                            # lower value
  z2 <- 0.9                            # higher value
  z3 <- 0.8
  z1<-z1^(1/b); z2<-z2^(1/b) ; z3 <- z3^(1/b)             # rescale variables

  rho <- FitnessGFEstimation(mc, death = death, plateff = plateff, model = model, muinf = muinf)

  if(rho$succeeds){

    if(init){
      mut <- MutationGFEstimation(mc, mfn = mfn, cvfn = cvfn, fitness = rho$fitness, death = death, plateff = plateff, model = model, muinf = muinf, init = init)
      c(mut, fitness = rho$fitness, succeeds = rho$succeeds)
    } else {
      a.est <- MutationGFEstimation(mc, fitness = rho$fitness, death = death, plateff = plateff, model = model, muinf = muinf, init = TRUE)
      Mutmodel <- new(FlanMutMod, list(
  		  mutations = a.est$mutations,
  		  fitness = rho$fitness,
  		  death = death,
  		  plateff = plateff,
  		  model = model,
        muinf = muinf,
        mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
  		))
      Cov <- Mutmodel$CovGFEstimation(z1, z2, z3)
      sd <- sqrt(Cov/length(mc))

      if(!is.null(mfn)) {
	if(cvfn > 0) {
	  if(plateff < 1) z3 <- 1-plateff+plateff*z3
	  pm.est <- Mutmodel$unbias.mutprob(sd[1], z3, mfn, cvfn)
	}
	else pm.est <- list(mutprob = a.est$mutations/mfn, sd.mutprob = sd[1]/mfn)

	c(pm.est, fitness = rho$fitness, sd.fitness = sd[2],
		  succeeds = rho$succeeds)

      } else list(mutations = a.est$mutations, sd.mutations = sd[1],
		  fitness = rho$fitness, sd.fitness = sd[2],
		  succeeds = rho$succeeds)
    }
  } else {
    res <- MutationGFEstimation(mc, mfn = mfn, cvfn = cvfn, fitness = rho$fitness, death = death, plateff = plateff, model = model, muinf = muinf, init = init)
    c(res, fitness = rho$fitness, sd.fitness = 0., succeeds = rho$succeeds)
  }

}

dclone <- function(m, fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }

  if(missing(model)){model <- "LD"}
  model <- match.arg(model)


  if(plateff < 0 | plateff > 1| length(plateff) > 1){
    stop("'plateff' must be a single positive and <= 1 number.")
  }

  if(model == "I"){
    if(!is.finite(muinf)){
      warning("'muinf' is infinite, 'model' is set to 'LD'.")
      model <- "LD"
    }
  } else {
    if(is.finite(muinf)) {
      warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      muinf <- Inf
    }
  }
  if(model == "H" & plateff < 1){
    # stop("If 'model' is 'H', 'plateff' can not be < 1.")
    warning("Probabilities are not available when 'model' is 'H' and 'plateff' < 1. 'plateff' is set to 1")
    plateff <- 1
  }

  # if(model == "LD") {
  #   integrands <- list(CLONE_P0_WD = function(x, rho, delta) {(1-x)*x^(rho-1)/(1-delta*x)},
	# 	     CLONE_PK_WD = function(x, rho, delta, k) {x^rho*(1-x)^(k-1)/(1-x*delta)^(k+1)}
	# 	     )
    # clone <- new(FlanExpClone, list(fitness = fitness, death = death, integrands = integrands))
  # }
  if(model == "LD") {
    # if(plateff < 1)
    clone <- new(FlanExpClone, list(fitness = fitness, death = death, plateff = plateff))
    # else clone <- new(FlanExpClone, list(fitness = fitness, death = death))
  } else if(model == "I") {
    # if(plateff < 1)
    clone <- new(FlanInhClone, list(fitness = fitness, death = death, plateff = plateff, muinf = muinf))
    # else clone <- new(FlanInhClone, list(fitness = fitness, death = death, muinf = muinf))
  } else clone <- new(FlanDirClone, list(fitness = fitness, death = death))


  M <- max(m)


  clone$dclone(M)[m+1]

}

dclone.dr <- function(m, fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }

  if(missing(model)){model <- "LD"}
  model <- match.arg(model)

  if(plateff < 0 | plateff > 1| length(plateff) > 1){
    stop("'plateff' must be a single positive and <= 1 number.")
  }

  if(model == "I"){
    if(!is.finite(muinf)){
      warning("'muinf' is infinite, 'model' is set to 'LD'.")
      model <- "LD"
    }
  } else {
    if(is.finite(muinf)) {
      warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      muinf <- Inf
    }
  }
  if(model == "H" & plateff < 1){
    # stop("If 'model' is 'H', 'plateff' can not be < 1.")
    warning("Probabilities are not available when 'model' is 'H' and 'plateff' < 1. 'plateff' is set to 1")
    plateff <- 1
  }
  # if(model == "LD") {
  #   integrands <- list(CLONE_P0_WD = function(x, rho, delta) {(1-x)*x^(rho-1)/(1-delta*x)},
	# 	     CLONE_PK_WD = function(x, rho, delta, k) {x^rho*(1-x)^(k-1)/(1-x*delta)^(k+1)}
	# 	     )
    # clone <- new(FlanExpClone, list(fitness = fitness, death = death, integrands = integrands))
  # }
  if(model == "LD") {
    # if(plateff < 1)
    clone <- new(FlanExpClone, list(fitness = fitness, death = death, plateff = plateff))
    # else clone <- new(FlanExpClone, list(fitness = fitness, death = death))
  } else if(model == "I") {
    # if(plateff < 1)
    clone <- new(FlanInhClone, list(fitness = fitness, death = death, plateff = plateff, muinf = muinf))
    # else clone <- new(FlanInhClone, list(fitness = fitness, death = death, muinf = muinf))
  } else clone <- new(FlanDirClone, list(fitness = fitness, death = death))


  M <- max(m)


  output <- clone$dclonedr(M)

  lapply(output, function(l) l[m+1])

}




dflan.grad <- function(m, mutations = 1., fitness = 1., death = 0., plateff = 1., model = c("LD", "H", "I"), muinf = +Inf, dalpha = TRUE, drho = FALSE){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }

  if(!(dalpha | drho)) stop("Need to precise which derivative has to be computed.")

  if(missing(model)){model = "LD"}
  model <- match.arg(model)

  if(plateff < 0 | plateff > 1| length(plateff) > 1){
    stop("'plateff' must be a single positive and <= 1 number.")
  }
  if(model == "I"){
    if(!is.finite(muinf)){
      warning("'muinf' is infinite, 'model' is set to 'LD'.")
      model <- "LD"
    }
    # if(plateff < 1) model <- "Ipef"
  } else {
    if(is.finite(muinf)) {
      warning(paste("if 'model' is '",model,"' 'muinf' is ignored.",sep=""))
      muinf <- Inf
    }
  }
  if(model == "H" & plateff < 1){
    # stop("If 'model' is 'H', 'plateff' can not be < 1.")
    warning("Probabilities are not available when 'model' is 'H' and 'plateff' < 1. 'plateff' is set to 1")
    plateff <- 1
  }
  # if(model == "LD" & plateff < 1) model <- "LDpef"

#   if(model == "LD") {
#     integrands <- list(CLONE_P0_WD = function(x, rho, delta) {(1-x)*x^(rho-1)/(1-delta*x)},
# 		       CLONE_PK_WD = function(x, rho, delta, k) {x^rho*(1-x)^(k-1)/(1-x*delta)^(k+1)},
# 		       CLONE_dP0_dr_WD = function(x, rho, delta) {(1-x)*x^(rho-1)/(1-delta*x)*log(x)},
# 		       CLONE_dPK_dr_WD = function(x, rho, delta, k) {x^rho*(1-x)^(k-1)/(1-x*delta)^(k+1)*log(x)}
# 		       )
#   } else integrands <- NULL
# #
  flan.mutmodel <- new(FlanMutMod, list(
    mutations = mutations,
    fitness = fitness,
    death = death,
    plateff = plateff,
    # integrands = integrands,
    model = model,
    muinf = muinf,
    mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
  ))

  M <- max(m)


  if(dalpha){
    if(drho) output <- flan.mutmodel$dflangrad(M)
    else output <- flan.mutmodel$dflanda(M)
  } else output <- flan.mutmodel$dflandr(M)

  lapply(output, function(l) l[m+1])
}


deduce.dflan <- function(m, mutations = 1., fitness = 1., death = 0., clone){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  M <- max(m)
#   if(length(clone) != M+1) {
#     stop("'clone' must be a vector with length 'max(m)'+1.")
#   }
#   if(missing(model)){model = "LD"}
#   model <- match.arg(model)


  flan.mutmodel <- new(FlanMutMod, list(
    mutations = mutations,
    fitness = fitness,
    death = death,
    model = "N", plateff = NULL, muinf = NULL, mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
  ))

  output <- flan.mutmodel$deduce.dflan(M, clone)

  output[m+1]
}

deduce.dflanda <- function(m, mutations = 1., fitness = 1., death = 0., clone){

  if(sum(m < 0) > 0 | sum(trunc(m) != m) > 0){
    stop("'m' must be a vector of positive integers")
  }
  if(mutations < 0 | length(mutations) > 1){
    stop("'mutations' must be a single positive number.")
  }
  if(fitness < 0 | length(fitness) > 1){
    stop("'fitness' must be a single positive number.")
  }
  if(death < 0 | death >= 0.5 | length(death) > 1){
    stop("'death' must be a single positive and < 0.5 number.")
  }
  M <- max(m)
#   if(length(clone) != M+1) {
#     stop("'clone' must be a vector with length 'max(m)'+1.")
#   }

#   if(missing(model)){model = "LD"}
#   model <- match.arg(model)

  flan.mutmodel <- new(FlanMutMod, list(
    mutations = mutations,
    fitness = fitness,
    death = death,
    model = "N", plateff = NULL, muinf = NULL, mc = NULL, mfn = NULL, cvfn = NULL, scale = NULL
  ))

  output <- flan.mutmodel$deduce.dflanda(M, clone)

  lapply(output, function(l) l[m+1])
}


## --------------------------------------------------------------------------------------------------------------
  ############################################## flantest object class ##############################################
  ## --------------------------------------------------------------------------------------------------------------

# Constructor of the class flantest
flantest <- function(Tstat, parameter, estimate, p.value, conf.int, estimates, null.value, alternative, model, method, data.name, sample.name, nsamples){
  obj <- list(Tstat = Tstat, estimate = estimate, parameter = parameter, conf.int = conf.int, p.value = p.value,
  null.value = null.value, alternative = alternative, data.name = data.name, model = model, method = method, nsamples = nsamples)
  class(obj) <- "flantest"
  obj
}

# Print function for flantest objects, inspired from print function for htest objects
print.flantest <- function(x, ...){
  digits <- getOption("digits")
  prefix <- "\t"
  cat("\n")
  if(x$nsamples == 2){
    cat(strwrap(paste("Two sample ", x$method, "-test", " (", x$model, " model)", sep = ""), prefix = prefix), sep = "\n")
  } else {
    cat(strwrap(paste("One sample ", x$method, "-test", " (", x$model, " model)", sep = ""), prefix = prefix), sep = "\n")
  }
  cat("\n")
  cat("--------------------------------- Data -------------------------------\n")

  if(x$nsamples == 2){
    if(length(x$data.name) == 1){
      cat("data:  ", x$data.name[[1]][1], " and ", x$data.name[[1]][2], "\n", sep = "")
    }
    if(length(x$data.name) == 2){
      if(length(x$data.name[[2]]) == 2){
	cat("data:  ", x$data.name[[1]][1], " with ", x$data.name[[2]][1],
	" and ", x$data.name[[1]][2], " with ", x$data.name[[2]][2], "\n", sep = "")
      } else {
	cat("data:  ", x$data.name[[1]][1], " with ", x$data.name[[2]][1],
	" and ", x$data.name[[1]][2], "\n", sep = "")
      }
    }
  }
  if(x$nsamples == 1){
    if(length(x$data.name) == 1){
      cat("data:  ", x$data.name[[1]][1], "\n", sep = "")
    }
    if(length(x$data.name) == 2){
      cat("data:  ", x$data.name[[1]][1], " with ", x$data.name[[2]][1], "\n", sep = "")
    }
  }

  out <- character()

  if(!is.null(x$null.value)){
    if(length(x$null.value) == 1){
      if(!is.null(x$parameter)){
	if(x$nsamples == 2){
	  Nparam <- dim(x$parameter)[2]
	  cat("Sample 1 parameters: ")
	  for(i in 1:Nparam){
      if(colnames(x$parameter)[i] == "mfn"){
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
        cat("                     ")
	      out <- c(paste(colnames(x$parameter)[i], " = ", format(x$parameter[1, i],
	      scientific = TRUE, digit = 5)))
	    } else {
	       out <- c(out, paste(colnames(x$parameter)[i], " = ", format(signif(x$parameter[1, i],
	      max(1, digits - 2)))))
      }
	  }
	  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
	  out <- character()
	  cat("Sample 2 parameters: ")
	  for(i in 1:Nparam){
      if(colnames(x$parameter)[i] == "mfn"){
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
        cat("                     ")
	      out <- c(paste(colnames(x$parameter)[i], " = ", format(x$parameter[2, i],
	      scientific = TRUE, digit = 5)))
	    } else {
	       out <- c(out, paste(colnames(x$parameter)[i], " = ", format(signif(x$parameter[2, i],
	      max(1, digits - 2)))))
      }
	  }
	} else {
	  Nparam <- length(x$parameter)
    cat("Sample parameters: ")
	  for(i in 1:Nparam){
	    if(names(x$parameter)[i] == "mfn"){
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
        cat("                   ")
	      out <- c(paste(names(x$parameter)[i], " = ", format(x$parameter[i],
	      scientific = TRUE, digit = 5)))
	    } else {
	      out <- c(out, paste(names(x$parameter)[i], " = ", format(signif(x$parameter[i],
	      max(1, digits - 2)))))
	    }
	  }
	}
      }
      cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
      cat("------------------------------ Statistics ----------------------------\n")
      if(!is.null(x$Tstat)){
	cat("Tstat = ", format(signif(x$Tstat,
	  max(1, digits - 2))), "\n")
      }
      if (!is.null(x$p.value)) {
	fp <- format.pval(x$p.value, digits = max(1, digits -
	    3))
	cat("p-value for", names(x$null.value[1]), if (substr(fp, 1, 1) ==
	  "<") fp else paste(" = ", fp), "\n")
      }

      if (!is.null(x$alternative)) {
	cat("Alternative hypothesis: ")
	alt.char <- switch(x$alternative, two.sided = "not equal to",
		less = "less than", greater = "greater than")
	cat("true ", names(x$null.value), " is ", alt.char,
	  " ", x$null.value, "\n", sep = "")
      }

      if (!is.null(x$conf.int)) {
	cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n",
	paste(format(x$conf.int), collapse = "  "),
	"\n")
      }
      if (!is.null(x$estimate)) {
	if(x$nsamples == 1){
	  cat("Sample estimates: \n")
	  print(x$estimate)
	}
	if(x$nsamples == 2){
	  cat("Sample 1 estimates: \n")
	  print(x$estimate[1, ])
	  cat("Sample 2 estimates: \n")
	  print(x$estimate[2, ])
	}
      }
    } else if(length(x$null.value) == 2){

      if(!is.null(x$parameter)){
	if(x$nsamples == 2){
	  cat("Sample 1 parameters: ")
	  Nparam <- dim(x$parameter)[2]
	  for(i in 1:Nparam){
      if(colnames(x$parameter)[i] == "mfn"){
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
        cat("                     ")
	      out <- c(paste(colnames(x$parameter)[i], " = ", format(x$parameter[1, i],
	      scientific = TRUE, digit = 5)))
	    } else {
	      out <- c(out, paste(colnames(x$parameter)[i], " = ", format(signif(x$parameter[1, i],
	      max(1, digits - 2)))))
	    }
	  }
	  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
	  out <- character()
	  cat("Sample 2 parameters: ")
	  for(i in 1:Nparam){
	    if(colnames(x$parameter)[i] == "mfn"){
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
        cat("                     ")
	      out <- c(paste(colnames(x$parameter)[i], " = ", format(x$parameter[2, i],
	      scientific = TRUE, digit = 5)))
	    } else {
	      out <- c(out, paste(colnames(x$parameter)[i], " = ", format(signif(x$parameter[2, i],
	      max(1, digits - 2)))))
	    }
	  }
	}
	if(x$nsamples == 1){
	  cat("Sample parameters: ")
	  Nparam <- length(x$parameter)
	  for(i in 1:Nparam){
      if(names(x$parameter)[i] == "mfn"){
        cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
        cat("                   ")
	      out <- c(paste(names(x$parameter)[i], " = ", format(x$parameter[i],
	      scientific = TRUE, digit = 5)))
	    } else {
	      out <- c(out, paste(names(x$parameter)[i], " = ", format(signif(x$parameter[i],
	      max(1, digits - 2)))))
	    }
	  }
	}
      }
      cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
      cat("------------------------------ Statistics ----------------------------\n")
      if(!is.null(x$Tstat)){
	cat("Tstat = (", format(signif(x$Tstat[1],
	  max(1, digits - 2))), " , ", format(signif(x$Tstat[2],
	  max(1, digits - 2))), ")", "\n", sep = "")
      }
      if (!is.null(x$p.value)) {
	fp1 <- format.pval(x$p.value[1], digits = max(1, digits -
	    3))
	fp2 <- format.pval(x$p.value[2], digits = max(1, digits -
	    3))
	cat("p-value for", names(x$null.value[1]), if (substr(fp1, 1, 1) ==
	  "<") fp1 else paste(" = ", fp1), "\np-value for", names(x$null.value[2]), if (substr(fp2, 1, 1) ==
	  "<") fp2 else paste(" = ", fp2), "\n")
      }
      if (!is.null(x$alternative)) {
	cat("Alternative hypotheses: ")
	alt.char1 <- switch(x$alternative[1], two.sided = "not equal to",
		less = "less than", greater = "greater than")
	alt.char2 <- switch(x$alternative[2], two.sided = "not equal to",
		less = "less than", greater = "greater than")
	cat("true ", names(x$null.value[1]), " is ", alt.char1,
	  " ", x$null.value[1], "\n", sep = "")
	cat("                        ")
	cat("true ", names(x$null.value[2]), " is ", alt.char2,
	  " ", x$null.value[2], "\n", sep = "")
      }
      if (!is.null(x$conf.int)) {
	cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval for ", names(x$null.value[1]), ": \n",
	  paste(format(x$conf.int[, 1]), collapse = "   "),
	  "\n",
	  format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval for ", names(x$null.value[2]), ": \n",
	  paste(format(x$conf.int[, 2]), collapse = "   "),
	  "\n", sep = "")
      }
      if (!is.null(x$estimate)) {
	if(x$nsamples == 1){
	  cat("Sample estimates: \n")
	  print(x$estimate)
	}
	if(x$nsamples == 2){
	  cat("Sample 1 estimates: \n")
	  print(x$estimate[1, ])
	  cat("Sample 2 estimates: \n")
	  print(x$estimate[2, ])
	}
      }
    }
    cat("\n")
    invisible(x)
  }
}


draw.clone <- function(t, mutprob = 1.e-2, fitness = 1., death = 0.,
	     dist = list("lnorm", meanlog = -0.3795851, sdlog = 0.3016223),
	     col = c("green4", "orange4")){
  #   Simulates a clone up to time t, represents the clone, and
  #   returns the vector of split times.
  #   At each division, the probability of having a mutant is mutprob.
  #   Variable death is the probability of death of a cell.
  #   The division times of have distribution dist, normalized
  #   to unit growth rate for mutant, to fitness for normal cells.
  #   The distribution dist is given as a list of a character chain
  #   followed by the list of parameters.
  #
  #   usage: draw.clone(5)
  #          draw.clone(t = 9, fitness = 0.5, death = 0.1)
  #          draw.clone(t = 3, mutprob = 0.1, fitness = 2, dist = dexp, death = 0.2)
  #

    if(!is.list(dist)){
      if(dist == "exp") {
        dist <- list("exp")
      } else if(dist == "dirac") {
        dist <- list("dirac")
      } else stop("'dist' must be a list of a character chain followed by its arguments if required")
    }
    names(dist)[1] <- "name"
    # if(dist$name == "exp"){
    #   names(dist)[2] <- "rate"
    # } else if(dist$name == "dirac"){
    #   names(dist)[2] <- "location"
    # } else
    if(dist$name == "lnorm"){
      names(dist)[2] <- "meanlog"
      names(dist)[3] <- "sdlog"
    } else if(dist$name == "gamma"){
      names(dist)[2] <- "shape"
      names(dist)[3] <- "scale"
    } else if(dist$name != "exp" & dist$name != "dirac") stop("'dist[[1]]' must be a character chain among 'exp', 'dirac', 'lnorm', 'gamma'")


		    # initialization
  p <- mutprob                            # mutation probability
  distn <- adjust.rate(dist, fitness = fitness, death = death)      # division times of normal cells
  distm <- adjust.rate(dist, fitness = 1., death = death)            # division times of mutants
  g <- 0                                  # current generation
  ce <- 1                                 # cells of current generation
  mu <- FALSE                                 # mutants in current generation
  bd <- 0                                 # birth dates of current generation
  dt <- rdt(1, distn)                      # life times of current generation
  de <- bd+dt                             # death dates of current generation
  or <- 1                                 # location of cells in current generation
  div <- {which((de<t)&                   # indices of cells dead before t
	      (runif(length(de))>death))} # that will divide
  de <- min(de, t)                         # truncate death dates
  nd <- length(div)                       # number of dividing cells
  nc <- nd                                # total number of cells
  CE <- ce                                # all cells
  MU <- mu                                # all mutants
  BD <- bd                                # all birth dates
  DE <- de                                # all death dates
  OR <- or                                # all locations
		    # main loop
  while((nd>0) & (nc<1e+4)){               # while divisions still happen
      g <- g+1                            # next generation
      ce <- 2*ce[div]; ce <- c(ce, ce+1)   # cells in next generation
	  mu <- mu[div]                   # mutants that divide
	  mu1 <- {sample(c(TRUE, FALSE), length(div),
	      replace = TRUE, prob = c(p, 1-p))}# new mutants in next generation
  # 	if (){mu1}
	  mu1 <- mu|mu1                   # mutants in next generation
	  mu <- c(mu1, mu)                 # all mutants in next generation
      nd <- nd*2                          # double the number of cells
      bd <- c(de[div], de[div])            # birth dates of daughters
	  imu <- which(mu)                # indices of mutants
	  nmu <- length(imu)              # number of mutants
	  ino <- which(!mu)               # indices of normal cells
	  nno <- length(ino)              # number of normal cells
      dtmu <- rdt(nmu, distm)              # division times of mutants
      dtno <- rdt(nno, distn)              # division times of normal cells
	  dt <- rep(0, nd)                 # division times of daughters
      dt[imu] <- dtmu                     # insert division times of mutants
	  dt[ino] <- dtno                 # insert division times of normal cells
      de <- bd+dt                         # death dates of daughters
      or<-c(or[div]-2^(-g), or[div]+2^(-g))# location of daughters
      div <- {which((de<t)&               # indices of cells dead before t
	      (runif(length(de))>death))} # that will divide
      de <- mapply(min, de, rep(t, nd))      # truncate death dates
      nd <- length(div)                   # number of dividing daughters
      nc <- nc+nd                         # total number of cells
      CE <- c(CE, ce)                      # stack new cells
      MU <- c(MU, mu)                      # stack mutants
      BD <- c(BD, bd)                      # stack birth dates
      DE <- c(DE, de)                      # stack death dates
      OR <- c(OR, or)                      # stack locations
  }                                       # end while
  if (nc > 1e+4){
    warning("too many cells to plot")
    }else{
		    # treatment for graphics

  colnor <- col[1]                      # color for normal cells
  colmut <- col[2]                   # color for mutants
  nc <- length(CE)                        # number of cells
  if (nc>1){                              # something to plot
  OR <- sort(OR, index.return = TRUE)$ix     # order for plotting
  CE <- CE[OR]                            # reorder cells
  MU <- MU[OR]                            # reorder mutants
  BD <- BD[OR]                            # reorder births
  DE <- DE[OR]                            # reorder deaths
  nmu <- length(which(MU))                # number of mutants
  nno <- nc - nmu                         # number of normal
  BDN <- BD[!MU]                          # births of normal
  DEN <- DE[!MU]                          # deaths of normal
  XN <- rbind(which(!MU), which(!MU))      # abscissas for vertical lines normal
  YN <- rbind(BDN, DEN)                    # ordinates for vertical lines normal
  nlc <- nc - length(which(DE<t))         # living cells at time t
  nlm <- nmu - length(which((DE<t)&MU))   # living mutants at time t
  para <- {substitute(list(               # list of parameters
	  "mutation probability" == t1, fitness == t2, death == t3,
	  cells == t4, mutants == t5),    # names
	  list(t1 = p, t2 = fitness, t3 = death,
	  t4 = nlc, t5 = nlm))}               # values
  name <- switch(dist[[1]],               # distribution of division times
	  dirac = "constant lifetimes",
	  exp = "exponential lifetimes",
	  gamma = "gamma lifetimes",
	  lnorm = "log-normal lifetimes")
  {matplot(XN, YN,                         # plot vertical lines
	  main = para,                     # add title
	  xlim = c(1, nc),
	  ylim = c(t, 0),                 # reverse y axis
	  type = "l", col = colnor,           # color for normal
	  lty = rep(1, nno),                # solid lines
	  axes = FALSE,                    # remove axes
	  xlab = name, ylab = "time")}        # axes labels
  axis(side = 2)                            # draw time axis
  if(nmu>0){                              # if any mutant
    BDM <- BD[MU]                        # births of mutant
    DEM <- DE[MU]                        # deaths of mutant
    XM <- rbind(which(MU), which(MU))     # abscissas for vertical lines mutant
    YM <- rbind(BDM, DEM)                 # ordinates for vertical lines mutant
    {matlines(XM, YM,                     # plot horizontal lines
	  type = "l", col = colmut,           # color for mutant
	  lty = rep(1, nmu))}               # solid lines
  }                                       # end if any mutants
  div <- which((DE<t)&(!MU))              # indices dividing cells normal
  DIN <- CE[div]                          # identities dividing cells normal
  DLN <- 2*DIN                            # identities of left daughters
  IO <- sort(CE, index.return = TRUE)        # numeral order of cells
  IOx <- IO$x; IOix <- IO$ix              # get cells and order
  indn <- which(IOx %in% DLN)             # identify left daughters normal
  XLN <- IOix[indn]                       # indices of left daughters normal
  XRN <- IOix[indn+1]                     # indices of right daughters normal
  XN <- rbind(XLN, XRN)                    # abscissas for horizontal lines
  YN <- rbind(BD[XLN], BD[XRN])            # ordinates for horizontal lines
  {matlines(XN, YN,                        # plot horizontal lines
	  type = "l", col = colnor,           # color normal
	  lty = rep(1, nno))}               # solid lines
  div <- which((DE<t)&MU)                 # indices dividing cells mutant
  if (length(div>0)){                     # if any dividing mutant
    DIM <- CE[div]                       # identities dividing cells mutant
    DLM <- 2*DIM                         # identities of left daughters
    indm <- which(IOx %in% DLM)          # identify left daughters mutant
    XLM <- IOix[indm]                    # indices of left daughters mutant
    XRM <- IOix[indm+1]                  # indices of right daughters mutant
    XM <- rbind(XLM, XRM)                 # abscissas for horizontal lines
    YM <- rbind(BD[XLM], BD[XRM])         # ordinates for horizontal lines
    {matlines(XM, YM,                     # plot horizontal lines
	  type = "l", col = colmut,           # color mutant
	  lty = rep(1, nmu))}               # solid lines
  }                                       # end if any dividing mutant
  }                                       # end if something to plot
  }                                       # end if graphics
}


rdt <- function(n, dist){
#   Returns a sample of size n of distribution dist.
#   The distribution dist is given as a list of a character chain
#   followed by the list of parameters.
#
#   usage: dist <- list("lnorm", meanlog = 0., sdlog = 1)
#          rdt(100, dist)
#          rdt(1000, goflnormBA)
#
  switch(dist[[1]],
      dirac = {                             # constant division times
	  l <- dist$location              # location parameter
	  rep(l, n)
	    },                            # end Dirac distribution
      exp = {                               # exponential distribution
	  r <- dist$rate                  # rate parameter
	  rexp(n, rate = r)
	    },                            # end exponential distribution
      gamma = {                             # gamma distribution
	  a <- dist$shape                 # shape parameter
	  l <- dist$scale                 # scale parameter
	  rgamma(n, shape = a, scale = l)
	    },                            # end gamma distribution
      lnorm = {                             # Log Normal distribution
	  a <- dist$sdlog                 # sdlog parameter
	  l <- dist$meanlog               # meanlog parameter
	  rlnorm(n, meanlog = l, sdlog = a)
	    }                             # end log normal distribution
  )                                       # end switch
}                                       # end function rdt


# rclone <- function(n, fitness = 1., death = 0., dist){
#
#   names(dist)[1] <- "name"
#   if(dist$name == "exp"){
#     names(dist)[2] <- "rate"
#   } else if(dist$name == "dirac"){
#     names(dist)[2] <- "location"
#   } else if(dist$name == "lnorm"){
#     names(dist)[2] <- "meanlog"
#     names(dist)[3] <- "sdlog"
#   } else if(dist$name == "gamma"){
#     names(dist)[2] <- "shape"
#     names(dist)[3] <- "scale"
#   } else stop("'dist[[1]]' must be a character chain among 'exp', 'dirac', 'lnorm', 'gamma'")
#
#
#   sim <- new(FlanSimClone, fitness, death, dist)
#
#   sim$rclone(n)
#
# }
#
# test.integral <- function(a = 0., b = 1., fitness = 1){
#
#   mInt <- new(FlanExpClone, list(fitness = fitness, death = 0))
#
#   mInt$integral(a, b)
# }
