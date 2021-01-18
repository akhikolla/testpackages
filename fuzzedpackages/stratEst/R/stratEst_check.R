#' Check model assumptions
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param model a fitted model of class \code{stratEst.model}.
#' @param chi.tests a logical. If \code{TRUE} chi square tests of global and local model fit are performed. Default is \code{FALSE}.
#' @param bs.samples an integer. The number of parametric bootstrap samples for the chi square tests. Default is 100.
#' @param verbose a logical, if \code{TRUE} messages of the checking process are printed to the console. Default is \code{FALSE}.
#' @export
#' @return A list of check results with the following elements:
#' \item{fit}{ a matrix. Contains the log likelihood, the number of free model parameters, and the value of the three information criteria.}
#' \item{chi.global}{a matrix. The results of the chi square test for global model fit.}
#' \item{chi.local}{ a matrix. The results of the chi square test for local model fit.}
#' @details The function for model checking of the package.
#' @references
#' Wang Z, Xu B, Zhou HJ (2014). "Social Cycling and Conditional Responses in the Rock-Paper-Scissors Game." \emph{Scientific Reports}, 4(1), 2045-2322.
#' @examples
#' ## Fit and check a mixture model for the rock-paper-scissors data of Wang, Xu, and Zhou (2014).
#' strategies.mixture = strategies.RPS[c("nash","imitate")]
#' model.mixture <- stratEst.model(data.WXZ2014,strategies.mixture)
#' model.mixture.check <- stratEst.check( model.mixture )
#' print(model.mixture.check$fit)
#' @export
stratEst.check <- function( model, chi.tests = F, bs.samples = 100, verbose = FALSE ){
  # check model
  if( "stratEst.model" %in% class(model) == F ){
    stop("stratEst.check error: The object passed to the argument 'model' must be of class 'stratEst.model'.")
  }

  # check test
  if ( "logical" %in% class(chi.tests) == F ){
    stop("stratEst.check error: The input argument 'chi.tests' must be a logical.");
  }

  # check bs.samples
  if ( bs.samples <= 0  | bs.samples%%1 != 0){
    stop("stratEst.check error: The number of bootstrap samples specified by the argument 'bs.samples' must be a positive integer. Default is 100.");
  }

  # check verbose
  if ( "logical" %in% class(verbose) == F ){
    stop("stratEst.check error: The input argument 'verbose' must be a logical.");
  }

  # retrieve info
  data <- model$fit.args$data
  shares <- model$shares
  LCR = is.null(model$coefficients) == F
  if( "list" %in% class(shares) ){
    num.strats <- length(shares[[1]])
  }else{
    num.strats <- length(shares)
  }
  num.ids = nrow(model$post.assignment)

  stratEst.check.return <- list()

  ###################################################################################
  # Global fit
  ###################################################################################

  fit <- matrix(model$fit[,1:5],1,5)
  colnames(fit) <- c("loglike","free.par","aic","bic","icl")
  rownames(fit) <- as.character(c(substitute(model)))
  stratEst.check.return$fit = fit

  ###################################################################################
  # Chi^2 tests
  ###################################################################################

  if( chi.tests ){
    # bootstrap the chi square test statistic
    data <- model$fit.args$data
    strategies <- model$strategies
    shares <- c(model$shares)
    sample.id <- model$fit.args$sample.id

    # data
    unique.choices <- unique(data$choice)

    greater <- rep(NA,bs.samples)
    chi <- rep(NA,bs.samples)
    greater_local <- matrix(NA,bs.samples,num.strats)
    chi_local <- matrix(NA,bs.samples,num.strats)
    if( verbose ){
      cat("boostrap chi^2 statistic\n")
    }

    for( m in 1:bs.samples ){
      sim.data <- stratEst.simulate(  data , strategies , shares , sample.id = sample.id )
      if( all(unique.choices %in% sim.data$choice) ){
        #selected strategies
        selected.strategies <- list()
        if( "list" %in% class(model$strategies[[1]]) ){
          est.strategies <- model$strategies[[1]]
        }else{
          est.strategies <- model$strategies
        }
        for( strs in 1:length(est.strategies) ){
          selected.strategies[[strs]] <- model$fit.args$strategies[[names(est.strategies)[strs]]]
        }
        names(selected.strategies) <- names(est.strategies)

        #selected shares
        selected.shares <- NULL
        if( is.null(model$fit.args$shares) == F ){
          if( "list" %in% class(model$fit.args$shares) ){
            selected.shares <- list()
            for( sam in 1:length(model$fit.args$shares) ){
              selected.shares.sample <- rep(NA,length(est.strategies))
              for( strs in 1:length(est.strategies) ){
                fit.shares.sample <- model$fit.args$shares[[sam]]
                selected.shares.sample[strs] <- fit.shares.sample[names(est.strategies)[strs] == names(model$fit.args$strategies)]
              }
              selected.shares[[sam]] <- selected.shares.sample
            }
            names(selected.shares) <- names(model$fit.args$shares)
          }else{
            selected.shares <- rep(NA,length(est.strategies))
            for( strs in 1:length(est.strategies) ){
              selected.shares[strs] <- model$fit.args$shares[names(est.strategies)[strs] == names(model$fit.args$strategies)]
            }
          }
        }

        # selected coefficients
        selected.coefficients <- NULL
        if( is.null(model$fit.args$coefficients) == F & length(model$strategies) > 1 ){
          selected.coefficients <- matrix(NA,nrow(model$fit.args$coefficients),length(model$strategies))
          for( strs in 1:length(model$strategies) ){
            selected.coefficients[,strs] <- model$fit.args$coefficients[,names(model$strategies)[strs]]
          }
          colnames(selected.coefficients) <- names(model$strategies)
          rownames(selected.coefficients) <- rownames(model$fit.args$coefficients)
        }

        selected.covariates <- NULL
        if( is.null(model$fit.args$covariates) == F & length(model$strategies) > 1 ){
          selected.covariates <- model$fit.args$covariates
        }

        # estimate model
        sim.model <- stratEst.model( sim.data , selected.strategies , selected.shares , coefficients = selected.coefficients, covariates = selected.covariates, sample.id = model$fit.args$sample.id, response = model$fit.args$response, sample.specific = model$fit.args$sample.specific, r.probs = model$fit.args$r.probs, r.trembles = model$fit.args$r.trembles, select = NULL, min.strategies = model$fit.args$min.strategies, crit = model$fit.args$crit, se = "analytic", outer.runs = model$fit.args$outer.runs, outer.tol = model$fit.args$outer.tol, outer.max = model$fit.args$outer.max, inner.runs = model$fit.args$inner.runs, inner.tol = model$fit.args$inner.tol, inner.max = model$fit.args$inner.max, lcr.runs = model$fit.args$lcr.runs, lcr.tol = model$fit.args$lcr.tol, lcr.max = model$fit.args$lcr.max, step.size = model$fit.args$step.size, penalty = model$fit.args$penalty, verbose = F )
        chi[m] <- sim.model$chi.global
        greater[m] <- sim.model$chi.global > model$chi.global
        chi_local[m,] <- sim.model$chi.local
        greater_local[m,] <- sim.model$chi.local > model$chi.local
      }
      if( verbose ){
        cat(paste("sample ",m," (of ",bs.samples,")\r",sep=""))
      }
    }
    if( verbose ){
     cat("\n")
    }

    p.value <- mean(greater,na.rm = T)
    p.values.local <- apply(greater_local,2,mean,na.rm = T)

    chi.global.result <- matrix(c(model$chi.global,min(chi,na.rm=T),mean(chi,na.rm=T),max(chi,na.rm=T),p.value),1,5)
    colnames(chi.global.result) <- c("chi^2","min","mean","max","p.value")
    rownames(chi.global.result) <- as.character(c(substitute(model)))
    stratEst.check.return$chi.global <- chi.global.result

    chi.local.result <- matrix(c(model$chi.local,apply(chi_local,2,min,na.rm=T),apply(chi_local,2,mean,na.rm=T),apply(chi_local,2,max,na.rm=T),p.values.local),num.strats,5)
    colnames(chi.local.result) <- c("chi^2","min","mean","max","p.value")
    rownames(chi.local.result) <- colnames(model$post.assignment)
    stratEst.check.return$chi.local <- chi.local.result

  }

  # ###################################################################################
  # # Non-Differential Measurement
  # ###################################################################################
  #
  # if( LCR & F ){
  #   measurement.result <- list()
  #   num.strats <- ncol(model$post.assignment)
  #   names.strategies <- colnames(model$post.assignment)
  #   id <- data$id
  #   choice <- data$choice
  #   choice.levels <- levels(choice)
  #   num.choice.levels <- length(choice.levels)
  #   num.ids <- nrow(model$post.assignment)
  #   unique.ids <- unique(id)
  #   covariates <- model$fit.args$covariates
  #   num.covariates <- length(covariates)
  #   covariate.mat <- matrix(NA,nrow(data),num.covariates)
  #   for( j in 1:num.covariates ){
  #     covariate.mat[,j] <- data[,covariates[j]]
  #   }
  #   num.strats <- ncol(model$post.assignment)
  #   hard.code <-  matrix(0,num.ids,num.strats)
  #   covariate <- matrix(NA,num.ids,num.covariates)
  #   colnames(covariate) <- covariates
  #   for( i in 1:num.ids ){
  #     max.indices <- which.max(model$post.assignment[ i , ])
  #     hard.code[ i , max.indices[1] ] = 1
  #     covariate[ i , ] <- unique( covariate.mat[data$id == unique.ids[i] , ] )
  #   }
  #   y <- as.factor(hard.code %*% c(1:num.strats))
  #   for( s in 1:num.strats ){
  #     ids.s <- unique.ids[ y == s ]
  #     choice.s <- choice[ id %in% ids.s ]
  #     covariates.s <- covariate.mat[ id %in% ids.s , ]
  #     if( length(unique(choice.s)) > 1 ){
  #       model.s.coefficients <- matrix(NA,num.choice.levels-1,num.covariates)
  #       colnames(model.s.coefficients) <- covariates
  #       rownames(model.s.coefficients) <- choice.levels[2:num.choice.levels]
  #       model.s.pvalue <- matrix(NA,num.choice.levels-1,num.covariates)
  #       colnames(model.s.pvalue) <- covariates
  #       rownames(model.s.pvalue) <- choice.levels[2:num.choice.levels]
  #       model.s <- nnet::multinom(choice.s ~ covariates.s - 1 , trace = F)
  #       model.s.coefficients <- matrix(summary(model.s)$coefficients,num.choice.levels-1,num.covariates)
  #       if( verbose ){
  #         cat(paste("boostrap ",names.strategies[s],"\n",sep=""))
  #       }
  #       for( m in 1:bs.samples ){
  #         ids.bs <- sample( ids.s , length(ids.s), replace = T )
  #         choice.bs <- NULL
  #         covariates.bs <- NULL
  #         for( i in 1:length(ids.bs) ){
  #           choice.bs <- c(choice.bs,choice[ id == ids.bs[i] ])
  #           covariates.bs <- rbind(covariates.bs,covariate.mat[ id == ids.bs[i] , ])
  #         }
  #         # FIX THIS!!!
  #         covariates.bs <- c(covariates.bs)
  #         if( length(unique(choice.bs)) > 1 ){
  #           model.bs <- nnet::multinom(choice.bs ~ covariates.bs - 1 , trace = F)
  #           coefs.bs <- summary(model.bs)$coefficients
  #           if( m == 1 ){
  #             model.s.pvalue = as.numeric( ( model.s.coefficients > 0 &  coefs.bs < model.s.coefficients ) | ( model.s.coefficients < 0 &  coefs.bs > model.s.coefficients ) )/bs.samples
  #           }else{
  #             model.s.pvalue = model.s.pvalue + as.numeric( ( model.s.coefficients > 0 &  coefs.bs < model.s.coefficients ) | ( model.s.coefficients < 0 &  coefs.bs > model.s.coefficients ) )/bs.samples
  #           }
  #         }
  #         if( verbose ){
  #           cat(paste("sample ",m," of (",bs.samples,")\r",sep=""))
  #         }
  #       }
  #       if( verbose ){
  #         cat("\n")
  #       }
  #       measurement.result[[s]] <- list( coefficients = model.s.coefficients , "p.values" = model.s.pvalue )
  #       if( verbose ){
  #         cat("\n")
  #       }
  #     }else{
  #       measurement.result[[s]] <- NA
  #     }
  #   }
  #   names(measurement.result) <- colnames(model$post.assignment)
  #   stratEst.check.return$measurement <- measurement.result
  # }
  #
  #
  #
  # ###################################################################################
  # # Regress posterior probability assignments on covariates
  # ###################################################################################
  #
  # if( LCR & F ){
  #   num.strats <- ncol(model$post.assignment)
  #   num.ids <- nrow(model$post.assignment)
  #   unique.ids <- unique(data$id)
  #   covariates <- model$fit.args$covariates
  #   num.covariates <- length(covariates)
  #   covariate.mat <- matrix(NA,nrow(data),num.covariates)
  #   for( j in 1:num.covariates ){
  #     covariate.mat[,j] <- data[,covariates[j]]
  #   }
  #   num.strats <- ncol(model$post.assignment)
  #   hard.code <-  matrix(0,num.ids,num.strats)
  #   covariate <- matrix(NA,num.ids,num.covariates)
  #   colnames(covariate) <- covariates
  #   for( i in 1:num.ids ){
  #     max.indices <- which.max(model$post.assignment[ i , ])
  #     hard.code[ i , max.indices[1] ] = 1
  #     covariate[ i , ] <- unique( covariate.mat[data$id == unique.ids[i] , ] )
  #   }
  #   y <- as.factor(hard.code %*% c(1:num.strats))
  #   levels(y) <- colnames(model$post.assignment)
  #   reference.level <- colnames(model$post.assignment)[1]
  #   y <- relevel(y, ref = reference.level)
  #   assignment.model <- nnet::multinom(y ~ covariate -1 , trace = F)
  #   assignment.model.summary <- summary(assignment.model)
  #   z <- abs(assignment.model.summary$coefficients/assignment.model.summary$standard.errors)
  #   p <- 2*stats::pt( z , model$res.degrees , lower = F )
  #   assignment.result = cbind( assignment.model.summary$coefficients , assignment.model.summary$standard.errors , z , rep(model$res.degrees,length(z)), p )
  #   colnames(assignment.result) <- c("coefficient","std.error","t-value","df","Pr(>|t|)")
  #   rownames(assignment.result) <- covariates
  #   stratEst.check.return$measurement <- assignment.result
  # }

  class(stratEst.check.return) <- c("stratEst.check","list")

  return(stratEst.check.return)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
