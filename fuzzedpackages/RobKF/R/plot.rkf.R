#' plot
#' @description A function to plot the output produced by \code{\link{AORKF_t}}, \code{\link{AORKF_huber}}, \code{\link{IORKF_huber}} or \code{\link{IOAORKF}}.
#' One can specify a time during the run for which the output should be displayed.
#' @name plot 
#' @param x An instance of an \code{ioaorkf} or \code{rkf} S3 class.
#' @param time A positive integer giving the time at which the output is to be displayed. It defaults to the number of observations.
#' @param subset A list of integers indicating the components of observations which are to be plotted.
#' @param conf_level A probability between 0 and 1 giving the confidence level at which the series are to be tested against anomalies. It defaults to 0.95.
#' @return A ggplot object.
#' @export
plot.rkf = function(x,time = NULL,subset = NULL,conf_level = 0.95,...){

  
   unexpectedarguments = names(list(...))
  
   if(length(unexpectedarguments)==1){warning(paste("The argument",unexpectedarguments,"has been ignored"))}
   if(length(unexpectedarguments)>1){warning(paste("The arguments",paste(unexpectedarguments,", "),"have been ignored"))}  
  
  is_observed<-value<-NULL
  
  if (is.null(time)){
    time = length(x[["Y"]])
  }
  
  time = as.integer(time)

  if (time > length(x[["Y"]])){
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
  
  mydaf= as.data.frame(t(Reduce(cbind,x[["Y"]])))
  n = nrow(mydaf)
  p = ncol(mydaf)
  
  if (is.null(subset)){
    subset = 1:p
  } 
  
  subset = unique(as.integer(subset))
  
  if (sum(subset %in% 1:p) < length(subset)){
    stop("subset has the wrong dimensions.")
  }
  
  colnames(mydaf) = paste("y",1:p,sep="")
  
  mydaf$x = 1:n
  
  mydaf = mydaf[,c("x",paste("y",subset,sep=""))]
  
  if (time < length(x[["Y"]])){
    
    molten_daf = molten.X<-melt(mydaf,id="x")
    
    molten.X$is_observed = 0
    molten.X[which(molten.X$x > time),"is_observed"] = 1
    names = as.vector(unique(molten.X$variable))
    
    out<-ggplot(data=molten.X)
    out<-out+aes(x=x,y=value,colour = is_observed)
    out<-out+theme(legend.position="none") 
    out<-out+ scale_colour_gradient(low="black", high="grey")
    out<-out+geom_point()
    out<-out+theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    
    out = out+geom_vline(xintercept = time+0.5)+ theme(legend.position="none")
    
    
  } else {
    
    molten_daf = molten.X<-melt(mydaf,id="x")
    names = as.vector(unique(molten.X$variable))
    
    out<-ggplot(data=molten.X)
    out<-out+aes(x=x,y=value)
    out<-out+theme(legend.position="none") 
    out<-out+geom_point()
    out<-out+theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    
  }
  
  if(length(pre_out)>0){
    
    for (ii in 1:length(pre_out)){
      
      if (x[["Type"]] == "IO"){
        
        out = out+geom_vline(xintercept = pre_out,colour="blue",alpha = 1)+ theme(legend.position="none")
        
      }
      
      if (x[["Type"]] == "AO"){
      
        for (jj in 1:length(subset)){
        
          Point_Anomalies = pre_out[ii] + (jj-1)*n
        
          out = out + geom_point(data = molten_daf[Point_Anomalies,] ,colour="red", size=1.5, alpha = 1)
        
        }
        
      }
      
    }
    
  }
  
  return(out)
  
}