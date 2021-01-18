#' Method dispatch for Generic Function Summary
#' @param object An object returned by the function\code{stratEst.check()}. An object of class \code{stratEst.check}.
#' @param ... additional arguments affecting the summary produced.
#' @export

summary.stratEst.check <- function( object , ...){

  if( "stratEst.check" %in% class(object) == F ){
    stop(paste("stratEst.summary error: The object ",as.character(object)," must be of class 'stratEst.check" ),sep="")
  }else{
    stratEst.check.return <- object
  }

  writeLines("==============================================================================================")
  writeLines("stratEst check summary")
  writeLines("==============================================================================================")
  writeLines("model fit:")
  print(stratEst.check.return$fit)
  writeLines("")
  if( is.null(stratEst.check.return$chi.global) == F ){
    writeLines("Chi^2 test of global model fit:")
    print(stratEst.check.return$chi.global)
    writeLines("")
  }
  if( is.null(stratEst.check.return$chi.local) == F ){
    writeLines("Chi^2 test of local model fit:")
    print(stratEst.check.return$chi.local)
    writeLines("")
  }
  writeLines("Please cite: Dvorak (2020). stratEst: strategy estimation in R.")
  writeLines("")


}
