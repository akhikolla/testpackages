.pass.class<-setClass("pass.class",representation(data="matrix",results="data.frame",Lmax="numeric",Lmin="numeric",alpha="numeric",lambda="numeric"))

pass.class<-function(data,results,Lmax,Lmin,alpha,lambda)
{
    .pass.class(data=data,results=results,Lmax=Lmax,Lmin=Lmin,alpha=alpha,lambda=lambda)
}

# not exported
pass.line.plot<-function(x,subset=1:ncol(x@data),variate_names=FALSE)
{
    # nulling out variables used in ggplot to get the package past CRAN checks
    k<-value<-NULL
    X<-as.data.frame(x@data[,subset])
    names<-paste("y",1:ncol(X),sep="")
    colnames(X)<-names
    # add in index variable for the time axis
    X<-cbind("k"=1:nrow(X),X)
    # melt the data
    molten.X<-melt(X,id="k")
    # generate multiple line plots ...
    p<-ggplot(data=molten.X)
    p<-p+aes(x=k,y=value)
    p<-p+geom_line()
    # check to see if there are any anomalies
    #if(!Reduce("||",is.na(x@results)))
    if (nrow(x@results) > 0)
    {
        p<-p+geom_vline(xintercept=x@results[,1],colour="red")
        p<-p+geom_vline(xintercept=x@results[,2],colour="red")
        p<-p+annotate("rect",xmin=x@results[,1],xmax=x@results[,2],ymin=-Inf,ymax=Inf,alpha=0.2,fill="yellow")
    }
    p<-p+facet_grid(variable~.,scales="free_y")
    # p<-p+facet_grid(factor(variable,levels=rev(names)) ~ .,scales="free_y")
    if(variate_names==FALSE)
    {
        p<-p+theme(axis.text.y=element_blank(),strip.text.y=element_blank())
    }
    return(p)
}

# not exported
pass.tile.plot<-function(x,subset=1:ncol(x@data),variate_names=FALSE)
{
    # nulling out variables used in ggplot to get the package past CRAN checks
    variable<-value<-NULL
    df<-as.data.frame(x@data[,rev(subset)])
    # normalise data
    for(i in 1:ncol(df))
    {
        df[,i]<-(df[,i]-min(df[,i]))/(max(df[,i])-min(df[,i]))
    }
    n<-data.frame("n"=seq(1,nrow(df)))
    molten.data<-melt(cbind(n,df),id="n")
    p<-ggplot(molten.data, aes(n,variable))
    p<-p+geom_tile(aes(fill=value))
    ymin<-0
    ymax<-ncol(df)
    # check to see if there are any anomalies
    #if(!Reduce("||",is.na(x@results)))
    if (nrow(x@results) > 0)
        {
            p<-p+annotate("rect",xmin=x@results[,1],xmax=x@results[,2],ymin=ymin,ymax=ymax+1,alpha=0.0,color="red",fill="yellow")
        }
    if(variate_names==FALSE)
    {
        p<-p+theme(axis.text.y=element_blank(),axis.title=element_blank())
    }
    return(p)
}

#' @name plot-pass.class
#'
#' @docType methods
#'
#' @rdname plot-methods
#'
#' @aliases plot,pass.class-method
#'
#' @export
setMethod("plot",signature=list("pass.class"),function(x,subset,variate_names=FALSE,tile_plot)
{
    if(missing(subset))
    {
        subset<-1:ncol(x@data)
    }
    if(missing(tile_plot))
    {
        tile_plot<-NULL
    }
    if(!is.logical(tile_plot))
    {
        if(is.null(tile_plot))
        {
            tile_plot<-FALSE
            if(ncol(as.matrix(x@data[,subset])) > 20)
            {
                tile_plot<-TRUE
            }
        }
        else
        {
            stop("tile_plot must be of type logical or NULL")
        }
    }
    if(!is.logical(variate_names))
    {
            stop("variable_names must be of type logical or NULL")
    }
    if(tile_plot)
    {
        p<-pass.tile.plot(x,subset,variate_names=FALSE)
    }
    else
    {
        p<-pass.line.plot(x,subset,variate_names)
    }
    return(p)
})

#' @name summary
#'
#' @docType methods
#' 
#' @rdname summary-methods
#'
#' @aliases summary,pass.class-method
#'
#' @export
setMethod("summary",signature=list("pass.class"),function(object,...)
{
    cat("PASS detecting change in mean","\n")
    cat("observations = ",nrow(object@data),"\n",sep="")
    cat("variates = ",ncol(object@data),"\n",sep="")
    cat("minimum segment length = ",object@Lmin,"\n",sep="")
    cat("maximum segment length = ",object@Lmax,"\n",sep="")
    cat("alpha = ",object@alpha,"\n",sep="") # tuning parameter
    cat("lambda = ",object@lambda,"\n",sep="")
    cat("Collective anomalies detected : ",nrow(object@results),"\n")
    invisible()
})


#' @name show
#'
#' @docType methods
#'
#' @aliases show,pass.class-method
#'
#' @rdname show-methods
#'
#' @export
setMethod("show",signature=list("pass.class"),function(object)
{
  summary(object)
  if (nrow(object@results) > 0){
    print(object@results)
  }
  cat("\n")
  invisible()
})


#' @name collective_anomalies
#'
#' @docType methods
#'
#' @rdname collective_anomalies-methods
#'
#' @aliases collective_anomalies,pass.class-method
#'
#' 
#' 
#' @export
setMethod("collective_anomalies",signature=list("pass.class"),function(object)
{
    return(object@results)
})
