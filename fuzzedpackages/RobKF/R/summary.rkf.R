
#' Summary
#' @name summary 
#' @description A function to summarise the output produced by \code{\link{AORKF_t}}, \code{\link{AORKF_huber}}, \code{\link{IORKF_huber}}, or \code{\link{IOAORKF}}.
#' One can specify a time during the run for which the output should be displayed.
#' @param object An instance of an \code{ioaorkf} or \code{rkf} S3 class.
#' @param time A positive integer giving the time at which the output is to be displayed. It defaults to the number of observations.
#' @param conf_level A probability between 0 and 1 giving the confidence level at which the series are to be tested against anomalies. It defaults to 0.95.
#' @export
summary.rkf = function(object,time = NULL,conf_level = 0.95,...){
  
  x = object
  
  unexpectedarguments = names(list(...))
  
  if(length(unexpectedarguments)==1){warning(paste("The argument",unexpectedarguments,"has been ignored"))}
  if(length(unexpectedarguments)>1){warning(paste("The arguments",paste(unexpectedarguments,", "),"have been ignored"))}
  
  if (is.null(time)){
    time = length(x[["Y"]])
  }
  
  time = as.integer(time)
  
  if (time > length(x[["Y"]]) ){
    stop("Time must be less than the number of observations.")
  }
  
  if (time < 1){
    stop("Time must be positive.")
  }
  
  conf_level = as.numeric(conf_level)
  
  if (conf_level >= 1){
    stop("conf_level must be between 0 and 1")
  }
  
  if (conf_level <= 0){
    stop("conf_level must be between 0 and 1")
  }
  
  scores = abs(Extract_all_anomalies(x))
  
  pre_out = which(scores> qchisq(conf_level,df = length(x[["Y"]][[1]])))
  
  pre_out = pre_out[which(pre_out <= time)]
  
  if (length(pre_out) == 0){
    
    cat(paste("At time",time,"out of", length(x[["Y"]]) , "no anomalies have been detected"))
    
  } else {
    
    cat(paste("At time",time,"out of", length(x[["Y"]]) , "the anomalies have been inferred at the following times:"))
    
    cat("\n")
    
    for (ii in 1:length(pre_out)){
      
      cat(paste(pre_out[ii] , "\t"))
      
    }
    
  }
  
}