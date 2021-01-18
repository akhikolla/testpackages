#' print
#'
#' @name print 
#' @param horizon A positive integer giving the smoothing horizon that is to be used. It must be at least equal to the number of rows of the horizonmatrix used to obtain the ioaorkf object.
#' @param ... Ignored.
#' @export
print.ioaorkf = function(x,time = NULL,horizon = NULL,...){
  
  unexpectedarguments = names(list(...))
  
  if(length(unexpectedarguments)==1){warning(paste("The argument",unexpectedarguments,"has been ignored"))}
  if(length(unexpectedarguments)>1){warning(paste("The arguments",paste(unexpectedarguments,", "),"have been ignored"))}
  
  if (is.null(horizon)){
    horizon = x$horizon
  }
  
  horizon = as.integer(horizon)
  
  if (horizon < x$horizon - 1){
    stop("horizon must be at least the number of rows of the horizonmatrix used to generate x.")
  }
  
  if (is.null(time)){
    time = length(x[["particles"]]) - 1
  }
  
  time = as.integer(time)
  
  if (time > length(x[["particles"]]) - 1){
    stop("Time must be less than the number of observations.")
  }
  
  if (time < 1){
    stop("Time must be positive.")
  }
  
  horizon = min(horizon,time)
  
  if (horizon > 0){
    x_new = Anomaly_Smoother(x,time,horizon)
  } else{ 
    x_new = x
  }
  
  pre_out = Extractanomalies(x_new)
  
  if (length(pre_out$anomaly_locations) == 0){
    
    cat(paste("At time",time,"out of", length(x_new[["particles"]]) - 1, "no anomalies have been detected"))
    
  } else {
    
    cat(paste("At time",time,"out of", length(x_new[["particles"]]) - 1, "the following anomalies have been inferred"))
    
    cat("\n")
    cat("\n")
    cat("location\t")
    cat("prob.\t\t")
    
    
    for (jj in 1:ncol(pre_out$anomaly_add_prob)){
      cat(paste("prob. add.", jj , "\t"))
    }
    
    for (jj in 1:ncol(pre_out$anomaly_inn_prob)){
      cat(paste("prob. inn.", jj , "\t"))
    }
    
    cat("\n")
    
    for (ii in 1:length(pre_out$anomaly_locations)){
      
      cat(paste(pre_out$anomaly_locations[ii] , "\t\t"))
      
      cat(paste(round(pre_out$anomaly_probs[ii],2) , "\t\t" ))
      
      
      for (jj in 1:ncol(pre_out$anomaly_add_prob)){
        cat(paste(round(pre_out$anomaly_add_prob[ii,jj],2) , "\t\t"))
      }
      
      for (jj in 1:ncol(pre_out$anomaly_inn_prob)){
        cat(paste(round(pre_out$anomaly_inn_prob[ii,jj],2) , "\t\t"))
      }
      
      cat("\n")
      
    }
  }
}