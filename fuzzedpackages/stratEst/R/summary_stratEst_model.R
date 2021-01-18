#' Method dispatch for Generic Function Summary
#' @param object An object returned by the estimation function\code{stratEst.model()}. An object of class \code{stratEst.model}.
#' @param ... additional arguments affecting the summary produced.
#' @export

summary.stratEst.model <- function( object , ... ){

  stratEst.return <- object

  convergence_string <- ifelse( is.null(stratEst.return$convergence) , "no parameters estimated" , ifelse( max( stratEst.return$convergence[ is.na(stratEst.return$convergence) == F ] ) < 0.001 , "yes" , ifelse( any( stratEst.return$convergence[ is.na(stratEst.return$convergence) == F ] < 0.001 ) , "partial" , "no") ) )

    writeLines("")
    writeLines(paste(deparse(substitute(object)),sep=""))
    writeLines(paste(rep("-",nchar(paste(deparse(substitute(object)),sep=""))),collapse = ""))
    writeLines(paste("number of individuals: ",stratEst.return$num.ids,sep=""))
    writeLines(paste("number of observations: ",stratEst.return$num.obs,sep=""))
    writeLines(paste("log likelihood: ",round(stratEst.return$loglike,2),sep=""))
    writeLines(paste("free model parameters: ",stratEst.return$free.par,sep=""))
    writeLines(paste("convergence: ",convergence_string,sep=""))
    writeLines("")
    writeLines("shares")
    writeLines(paste(rep("-",nchar("shares")),collapse = ""))
    if( "list" %in% class(stratEst.return$shares) ){
      print(round(do.call(rbind,stratEst.return$shares),2))
    }else{
      print(round(stratEst.return$shares,2))
    }
    writeLines("")
    if( is.null(stratEst.return$coefficients) == F ){
      writeLines("latent class coefficients")
      writeLines(paste(rep("-",nchar("latent class coefficients")),collapse = ""))
      print(round(stratEst.return$coefficients,2))
      writeLines("")
    }
    writeLines("strategies")
    writeLines(paste(rep("-",nchar("strategies")),collapse = ""))
    if( "list" %in% class(stratEst.return$strategies[[1]])  ){
      strategies_sample_list <- NULL
      strategies_print <- stratEst.return$strategies
      names_samples <- names(stratEst.return$strategies)
      for( i in 1:length(stratEst.return$strategies) ){
        strategies_sample <- do.call(rbind,strategies_print[[i]])
        row_names_strategies_sample <- rownames(strategies_sample)
        rownames(strategies_sample) =  paste( names_samples[i], "." , row_names_strategies_sample, sep="")
        strategies_sample_list <- rbind( strategies_sample_list , strategies_sample )
      }
      print(strategies_sample_list,2)
    }else{
      strategies_print <- do.call(rbind,stratEst.return$strategies)
      print(round(strategies_print,2))
    }
    writeLines("")
}




# summary.stratEst.model <- function( object , ... , objects = c("model","fit","shares","coefficients","strategies","parameters")){
#
#   stratEst.return <- object
#
#   if( "character" %in% class(objects) == F ){
#     stop(paste("stratEst.summary error: The object ",as.character(objects)," supplied as argument 'objects' must be a character vector.",sep=""))
#   }else{
#     for( i in 1:length(objects)){
#       if( objects[i] %in% c("model","fit","shares","coefficients","strategies","parameters") == F ){
#         stop(paste("stratEst.summary error: The object ",as.character(objects)," can only contain the character strings: model, fit, shares, coefficients, strategies, and parameters.",sep=""))
#       }
#     }
#   }
#   convergence_string <- ifelse( is.null(stratEst.return$convergence) , "no parameters estimated" , ifelse( max( stratEst.return$convergence[ is.na(stratEst.return$convergence) == F ] ) < 0.001 , "yes" , ifelse( any( stratEst.return$convergence[ is.na(stratEst.return$convergence) == F ] < 0.001 ) , "partial" , "no") ) )
#
#   if( all( objects == c("model","fit","shares","coefficients","strategies","parameters") ) ){
#     writeLines("==============================================================================================")
#     writeLines("stratEst summary")
#     writeLines("==============================================================================================")
#   }
#   if( "model" %in% objects ){
#     writeLines("model:")
#     writeLines(paste(rep("-",nchar("model:")),collapse = ""))
#     writeLines(paste("number of individuals: ",stratEst.return$num.ids,sep=""))
#     writeLines(paste("number of observations: ",stratEst.return$num.obs,sep=""))
#     writeLines(paste("number of model parameters: ",stratEst.return$num.par,sep=""))
#     writeLines(paste("number of free model parameters: ",stratEst.return$free.par,sep=""))
#     writeLines(paste("residual degrees of freedom: ", stratEst.return$res.degrees,sep=""))
#     writeLines(paste("convergence: ",convergence_string,sep=""))
#     writeLines("")
#   }
#   if( "fit" %in% objects ){
#     writeLines("fit:")
#     writeLines(paste(rep("-",nchar("fit:")),collapse = ""))
#     writeLines(paste("log likelihood: ",round(stratEst.return$loglike,2),sep=""))
#     writeLines(paste("model entropy: ",round(stratEst.return$entropy.model,2),sep=""))
#     writeLines(paste("aic: ",round(stratEst.return$aic,2),sep=""))
#     writeLines(paste("bic: ",round(stratEst.return$bic,2),sep=""))
#     writeLines(paste("icl: ",round(stratEst.return$icl,2),sep=""))
#     writeLines("")
#   }
#   if( "shares" %in% objects ){
#     writeLines("shares:")
#     writeLines(paste(rep("-",nchar("shares:")),collapse = ""))
#     if( "list" %in% class(stratEst.return$shares) ){
#       #writeLines(paste(rep("=",(sum(nchar(colnames(stratEst.return$shares[[1]])))+max(nchar(rownames(stratEst.return$shares[[1]])))+length(stratEst.return$shares[[1]]))), collapse = ""))
#       print(round(do.call(rbind,stratEst.return$shares),2))
#     }else{
#       print(stratEst.return$shares,2)
#     }
#     writeLines("")
#   }
#   if( "coefficients" %in% objects ){
#     if( is.null(stratEst.return$coefficients) == F ){
#       writeLines("latent class coefficients:")
#       writeLines(paste(rep("-",nchar("latent class coefficients:")),collapse = ""))
#       print(round(stratEst.return$coefficients,2))
#       writeLines("")
#     }
#   }
#   if( "strategies" %in% objects ){
#     writeLines("strategies:")
#     writeLines(paste(rep("-",nchar("strategies:")),collapse = ""))
#     if( "list" %in% class(stratEst.return$strategies[[1]])  ){
#       strategies_sample_list <- NULL
#       strategies_print <- stratEst.return$strategies
#       names_samples <- names(stratEst.return$strategies)
#       for( i in 1:length(stratEst.return$strategies) ){
#         strategies_sample <- do.call(rbind,strategies_print[[i]])
#         row_names_strategies_sample <- rownames(strategies_sample)
#         rownames(strategies_sample) =  paste( names_samples[i], "." , row_names_strategies_sample, sep="")
#         strategies_sample_list <- rbind( strategies_sample_list , strategies_sample )
#       }
#       print(strategies_sample_list,2)
#     }else{
#       strategies_print <- do.call(rbind,stratEst.return$strategies)
#       print(round(strategies_print,2))
#     }
#     writeLines("")
#   }
#   # if( "parameters" %in% objects ){
#   #   writeLines("parameters:")
#   #   writeLines(paste(rep("-",nchar("parameters:")),collapse = ""))
#   #
#   #   par_matrix <- NULL
#   #
#   #   if( length(stratEst.return$shares.par) > 1 & is.null(stratEst.return$coefficients)  ){
#   #     par <- stratEst.return$shares.par
#   #     se <- stratEst.return$shares.se
#   #     se[ se == 0 ] = NA
#   #     z <- abs(par/se)
#   #     p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
#   #     share_matrix = cbind( par , stratEst.return$shares.quantiles , se , z , p )
#   #     colnames(share_matrix) <- c("estimate",colnames(stratEst.return$shares.quantiles),"std.error","t-value","Pr(>|t|)")
#   #     rownames(share_matrix) <- paste("shares.par.",as.character(seq(1,nrow(share_matrix),by = 1)),sep="")
#   #     par_matrix = rbind(par_matrix,share_matrix)
#   #   }
#   #   if( is.null(stratEst.return$coefficients.par) == F ){
#   #     par <- stratEst.return$coefficients.par
#   #     se <- stratEst.return$coefficients.se
#   #     se[ se == 0 ] = NA
#   #     z <- abs(par/se)
#   #     p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
#   #     coefficients_matrix = cbind( par ,  stratEst.return$coefficients.quantiles , se , z , p )
#   #     colnames(coefficients_matrix) <- c("estimate",colnames(stratEst.return$coefficients.quantiles),"std. error","t-value","Pr(>|t|)")
#   #     rownames(coefficients_matrix) <- paste("coefficients.par.",as.character(seq(1,nrow(coefficients_matrix),by = 1)),sep="")
#   #     par_matrix = rbind(par_matrix,coefficients_matrix)
#   #
#   #   }
#   #   if( length(stratEst.return$probs.par > 0) & length(stratEst.return$probs.se > 0) ){
#   #     par <- stratEst.return$probs.par
#   #     se <- stratEst.return$probs.se
#   #     se[ se == 0 ] = NA
#   #     z <- abs(par/se)
#   #     p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
#   #     response_matrix = cbind( par , stratEst.return$probs.quantiles , se , z , p )
#   #     colnames(response_matrix) <- c("estimate",colnames(stratEst.return$probs.quantiles),"std.error","t-value","Pr(>|t|)")
#   #     rownames(response_matrix) <- paste("probs.par.",as.character(seq(1,nrow(response_matrix),by = 1)),sep="")
#   #     par_matrix = rbind(par_matrix,response_matrix)
#   #   }
#   #   if( length(stratEst.return$trembles.par > 0) ){
#   #     par <- stratEst.return$trembles.par
#   #     se <- stratEst.return$trembles.se
#   #     se[ se == 0 ] = NA
#   #     z <- abs(par/se)
#   #     p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
#   #     tremble_matrix = cbind( par , stratEst.return$trembles.quantiles , se , z , p )
#   #     colnames(tremble_matrix) <- c("estimate",colnames(stratEst.return$trembles.quantiles),"std.error","t-value","Pr(>|t|)")
#   #     rownames(tremble_matrix) <- paste("trembles.par.",as.character(seq(1,nrow(tremble_matrix),by = 1)),sep="")
#   #     par_matrix = rbind(par_matrix,tremble_matrix)
#   #   }
#   #
#   #   if( is.null(par_matrix) == F ){
#   #     print(round(par_matrix,3))
#   #     writeLines("")
#   #   }
#   #
#   # }
#   # if( all( objects == c("model","fit","shares","coefficients","strategies","parameters") ) ){
#   #   writeLines("Please cite: Dvorak (2020). stratEst: strategy estimation in R.")
#   #   writeLines("")
#   # }
#
#
# }
