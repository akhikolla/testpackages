


.changepoint.mv.class<-setClass("changepoint.mv.class",representation(data="data.frame",cpts.uv="list",cpts.mv="list",alpha="numeric",cost="character",mad="logical"))

changepoint.mv.class<-function(data,cpts.uv,cpts.mv,alpha,cost,mad,...)
{
    .changepoint.mv.class(data=data,cpts.uv=cpts.uv,cpts.mv=cpts.mv,alpha=alpha,cost=cost,mad=mad,...)
}


.changepoint.mv.mrc.class<-setClass("changepoint.mv.mrc.class",contains="changepoint.mv.class",representation(pmax="numeric",penalised.costs="numeric",locations="list",affected="list"))

changepoint.mv.mrc.class<-function(data,cpts.uv,cpts.mv,alpha,cost,mad,min.dist,pmax,penalised.costs,locations,affected,...)
{
    .changepoint.mv.mrc.class(changepoint.mv.class(data=data,cpts.uv=cpts.uv,cpts.mv=cpts.mv,alpha=alpha,cost=cost,mad=mad),pmax=pmax,
                              penalised.costs=penalised.costs,locations=locations,affected=affected,...)
}





if(!isGeneric("cpts.uv")) {setGeneric("cpts.uv",function(x,...) {standardGeneric("cpts.uv")})}
#' Univariate changepoint locations.
#'
#' @name cpts.uv
#'
#' @description Returns a list of vectors containing the univariate changepoint locations.
#'
#' @param x An S4 object as returned by \code{\link{mrc}}.
#'
#' @return A list of \eqn{N} vectors containing the univariate changepoint locations. Each vector corresponds to an individual variate in the data.
#'
#' @docType methods
#'
#' @aliases cpts.uv,changepoint.mv.mrc.class-method
#'
#' @rdname cpts.uv-methods
#'
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample)
#' cpts.uv(res)
#' 
#' @export
setMethod("cpts.uv",signature=list("changepoint.mv.mrc.class"),
          function(x)
          {
              if(length(x@cpts.uv)==0)
              {
                  stop("univariate results unavailable\n")
              }
              x@cpts.uv
          }
          )


if(!isGeneric("cpts.mv")) {setGeneric("cpts.mv",function(x,...) {standardGeneric("cpts.mv")})}
#' Multivariate changepoint locations.
#'
#' @name cpts.mv
#'
#' @description Returns a list of vectors containing the multivariate changepoint locations.
#'
#' @param x An S4 object as returned by \code{\link{mrc}}.
#' @param p The number of most recent changepoints locations to be consisdered. Default value is \code{p=x@pmax} where pmax
#' 
#' @return A list of \eqn{N} vectors containing the multivariate changepoint locations. Each vector corresponds to an individual variate in the data.
#'
#' @docType methods
#'
#' @aliases cpts.mv,changepoint.mv.mrc.class-method
#'
#' @rdname cpts.mv-methods
#'
#'
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample)
#' cpts.mv(res)
#' cpts.mv(res,p=3)
#'
#' @export
setMethod("cpts.mv",signature=list("changepoint.mv.mrc.class"),
          function(x,p=NULL)
          {
              if(is.null(p))
              {
                  p<-MDL(x)
              }
              assert_is_numeric(p)
              if(p>0 && p<=x@pmax)
                  {
                      return(x@cpts.mv[[p]])
                  }
              cat("p must be a positive integer with 0 < p <= pmax\n")
          }
          )



if(!isGeneric("cpts.mr")) {setGeneric("cpts.mr",function(x,...) {standardGeneric("cpts.mr")})}
#' Most recent changepoint locations.
#'
#' @name cpts.mr
#'
#' @description Returns a list of vectors containing the most recent changepoint locations.
#'
#' @param x An S4 object as returned by \code{\link{mrc}}.
#' @param p The number of most recent changepoints locations to be considered. Default value is \code{p=x@pmax} where pmax
#' is the value specified when \code{\link{mrc}} was called.
#'
#' @return A data frame containing the most recent changepoint locations and the the variates corresponding to those locations.
#'
#' @docType methods
#'
#' @aliases cpts.mr,changepoint.mv.mrc.class-method
#'
#' @rdname cpts.mr-methods
#'
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample)
#' cpts.mr(res)
#' cpts.mr(res,p=2)
#'
#' @export
setMethod("cpts.mr",signature=list("changepoint.mv.mrc.class"),
          function(x,p=x@pmax)
          {
	      if(p <= 0 || p > x@pmax)
	      {
	         stop("p must be in the range 0 < p <= pmax")
	      }
              mrcpts<-data.frame("cpt"=unlist(Map(function(var) var[length(var)-1],cpts.mv(x,p))),"variable"=1:(ncol(x@data)-1))
              mrcpts<-Reduce(rbind,split(mrcpts,mrcpts$cpt))
              row.names(mrcpts)<-1:nrow(mrcpts)
              return(mrcpts)
          }
          )



if(!isGeneric("data.set")) {setGeneric("data.set",function(x,...) {standardGeneric("data.set")})}
#' Recovers the data from the results of a changepoint analysis.
#'
#' @name data.set
#'
#' @description Recovers the data from the S4 class returned by  \code{\link{mrc}}. The data is stored in a data frame and, unless \code{indexed=TRUE}
#' was specified when the data was analysed, the first column will contain an index variable. If the original data did not have column names, default ones of the form V.n, where n is the column
#' number, will be added. 
#'
#' @param x An S4 class instance obtained from \code{\link{mrc}}.
#'
#' @return A data frame containing the original data and an index variable (if \code{indexed=TRUE} was used for the analysis).
#'
#' @docType methods
#'
#' @aliases data.set,changepoint.mv.mrc.class-method
#'
#' @rdname data.set-method
#'
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample)
#' head(data.set(res))
#' 
#' @export
setMethod("data.set",signature=list("changepoint.mv.mrc.class"),
          function(x)
          {
              x@data
          }
          )




#' Displays S4 objects produced by changepoint.mv methods
#'
#' @name show
#'
#' @description Displays S4 object produced by \code{\link{mrc}}. The information produced is the same as that provided by the \code{summary} method.
#' The method is used by the S4 system for automatic printing.
#'
#' @param object An S4 object produced by \code{\link{mrc}}.
#'
#' @docType methods
#'
#' @aliases show,changepoint.mv.mrc.class-method
#'
#' @rdname show-methods
#'
#'@examples
#' 
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample)
#' # the following lines all produce the same output
#' res
#' summary(res)
#' show(res)
#' print(res)
#'
#' @export
setMethod("show",signature=list("changepoint.mv.mrc.class"),function(object)
{
    summary(object)
})

#' Summary
#'
#' @name summary
#'
#' @description Produces and prints a summary of the results contained in an S4 class instance produced by \code{\link{mrc}}. The information
#' produced depends on the analysis type and options used, but typically includes the total penalised cost, penalty values and cost function type.
#'
#' @param object An S4 class instance produced by \code{\link{mrc}}.
#'
#' @docType methods
#'
#' @rdname summary-method
#'
#' @aliases summary,changepoint.mv.class-method
#'
#' @export
setMethod("summary",signature=list("changepoint.mv.class"),function(object)
{
    
    if(length(object@cpts.uv)>0 && length(object@cpts.mv) > 0)
    {
        cat("Univariate and Multivariate Analysis","\n")
    }
    if(length(object@cpts.uv)>0 && length(object@cpts.mv) == 0)
    {
        cat("Univariate Analysis","\n")
    }
    if(length(object@cpts.uv) == 0 && length(object@cpts.mv) > 0)
    {
        cat("Multivariate Analysis","\n")
    }
    #if(length(object@cpts.uv)>0)
    #{
    #    cat("Univariate Changepoints :","\n")
    #    Map(function(varname,cpts) cat(varname," : ",cpts,"\n"),names(object@data)[-1],object@cpts.uv)
    #}
    #if(length(object@cpts.mv)>0)
    #{
    #    cat("Multivariate Changepoints :","\n")
    #    Map(function(varname,cpts) cat(varname," : ",cpts,"\n"),names(object@data)[-1],object@cpts.mv)
    #}
    cat("Cost Function : ",object@cost,"\n")
    # cat("Penalty : ",object@penalty,"\n")
    cat("Minimum Distance between Changepoints : ",object@min.dist,"\n")
    cat("Alpha Value : ",object@alpha,"\n")    
    invisible()
})


if(!isGeneric("costs")) {setGeneric("costs",function(x,...) {standardGeneric("costs")})}
#' Uses the changepoint locations to determine the penalised cost of the segmented data.
#' 
#' @name costs
#'
#' @description
#' For results obtained using \code{\link{mrc}}, \code{costs} calculates the the total penalised cost for all segments across all variates for different numbers of
#' most recent change point values (p). It also calculates the code length \eqn{{\rm log}_{2}n^{p}p^{N}} and the total of the penalised cost and code length. The result
#' is a data frame containing p, the penalised cost, the code length and the sum of penalised cost and code length. The row in the data frame with the smallest total corresponds to
#' the minimum description length (MDL). See Bardwell, Eckley, Fearnhead and Smith, (2016) for more details about the cost and code length.
#'
#' @param x An S4 object as returned by \code{\link{mrc}}.
#' 
#' @return Data frame containing cost information as described in the description section (above).
#'
#' @docType methods
#'
#' @aliases costs,changepoint.mv.mrc.class-method
#'
#' @rdname costs-methods
#'
#' @references \insertRef{doi:10.1080/00401706.2018.1438926}{changepoint.mv}
#'
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample,pmax=8)
#' costs(res)
#'
#' @export
setMethod("costs",signature=list("changepoint.mv.mrc.class"),function(x)
{
    n<-nrow(x@data)
    N<-ncol(x@data)
    code.length<-unlist(Map(function(p) N*log(p,2) + p*log(n,2),1:x@pmax))
    return(as.matrix(data.frame("p"=1:length(x@penalised.costs),
                                "cost"=x@penalised.costs,
                                "code.length"=code.length,
                                "total"=x@penalised.costs+code.length)))
})





if(!isGeneric("MDL")) {setGeneric("MDL",function(x,...) {standardGeneric("MDL")})}
#' Calculates the Minimum Description Length.
#'
#' @name MDL
#'
#' @description Calculates the Minimum Description Length (MDL) using the result obtained from \code{\link{mrc}}. The MDL indicates how many most recent changepoints there are for the data.
#' For a full definition of the MDL and a description of its calculation see Bardwell, Eckley, Fearnhead and Smith, (2016). 
#'
#' @param x An S4 object as returned by \code{\link{mrc}}.
#'
#' @return The Minimum Description Length (MDL).
#' 
#' @docType methods
#'
#' @rdname MDL-methods
#'
#' @aliases MDL,changepoint.mv.mrc.class-method
#'
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample,pmax=8)
#' MDL(res)
#'
#' @seealso \code{\link{mrc}}
#'
#' @references \insertRef{doi:10.1080/00401706.2018.1438926}{changepoint.mv}
#'
#' @export
setMethod("MDL",signature=list("changepoint.mv.mrc.class"),function(x)
{
    return(which.min(x@penalised.costs + ((ncol(x@data)-1)*log(1:x@pmax)+ (1:x@pmax)*log(nrow(x@data)))/log(2))) ##MDL    
})



plot.heatmap<-function(x,multivariate=TRUE,variable.names=FALSE)
{
# nulling out variables used in ggplot to get the package past CRAN checks
k<-variable<-value<-cpt<-NULL # note, the value of NULL will not really be used - thsi is a hack to satisfy CRAN checks
# ggplot heatmap plots "upside down" - so reverse the columns of data 
df<-cbind(x@data[,1,drop=FALSE],rev(x@data[,2:ncol(x@data)]))
if(multivariate)
{
    df.cpts<-Reduce(cbind,Map(function(x) {res<-rep(0,nrow(df));res[x[-1]]<-1;data.frame(res);},x@cpts.mv))
}
else
{
    df.cpts<-Reduce(cbind,Map(function(x) {res<-rep(0,nrow(df));res[x[-1]]<-1;data.frame(res);},x@cpts.uv))
}
# need to reverse the changepoints as well 
df.cpts<-rev(df.cpts)
df.cpts<-cbind(df[,1],df.cpts)
df[,2:ncol(df)]<-apply(df[,2:ncol(df)],2,function(x) (x-min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))
df<-melt(df,id.var="k",value.name="value")
names(df.cpts)<-c("k",names(x@data)[-1])
df.cpts<-melt(df.cpts,id.var="k",value.name="cpt")
df<-cbind(df,df.cpts$cpt)
names(df)[4]<-"cpt"
p<-ggplot(df,aes(k,variable))
p<-p+geom_tile(aes(fill=value,colour=as.factor(cpt)),width=1,height=1,size=1)
p<-p+scale_colour_manual(values=alpha(c("black","yellow"),c(0,1)))
p<-p+guides(colour=FALSE)
if(variable.names==FALSE)
{
    p<-p+theme(axis.text.y=element_blank())
}

return(p)    

}



#' Visualisation of data and changepoint locations.
#'
#' @name plot
#'
#' @description Plot methods for S4 objects returned by \code{\link{mrc}}. The plot produced depends on the type of the
#' S4 object.
#'
#' For objects produced by \code{\link{mrc}} a heatmap of the data is displayed along with the location of the univariate changepoints (yellow), most recent univariate changepoints (green),
#' and most recent multivariate changepoint locations (red). 
#'
#' A number of arguments with default values are provided to control aspects of how the data and changepoint locations are displayed. The plot methods return a ggplot2 object which can be
#' modified if required.
#'
#' @docType methods
#'
#' @param x An S4 class returned by \code{\link{mrc}}. 
#' @param p Integer value indicating the number of most recent changepoint locations to consider. The minimum value is 1 and the maximum value is the value of \code{pmax} used
#' in the call to \code{\link{mrc}}. Default value is the \code{p=MDL(x)} (see \code{\link{MDL}} for further details).
#' @param group Logical value used to indicate if the variates that share a most recent changepoint should be grouped together in the plot. Default is \code{group=FALSE}.
#' @param display.variable.names Logical value. If \code{display.variable.names=TRUE} then the variable names are displayed in the plot. If there are a large number of variates in the data, then it
#' can be useful disable variable names by setting \code{display.variable.names=FALSE}. Default is \code{display.variable.names=TRUE}.
#' @param show Logical value used to indicate if the plot should be displayed (the ggplot object produced by plot is always invisibly returned). Default value is \code{show=TRUE}.
#' 
#' @return Invisibly returns a ggplot object.
#'
#' @rdname plot-methods
#'
#' @aliases plot,changepoint.mv.mrc.class,ANY-method
#'
#'
#'
#' @examples
#' \dontrun{
#' # visualising most recent changepoints 
#' data(mrcexample)
#' res<-mrc(mrcexample[,1:10])
#' p.1<-plot(res,p=2)
#' p.2<-plot(res,p=5)
#' p.3<-plot(res,p=2,group=TRUE)
#' p.4<-plot(res,p=5,group=TRUE)
#' if(require(gridExtra))
#' {
#'   grid.arrange(p.1,p.2,p.3,p.4)
#' }
#' }
#' @export
setMethod("plot",signature=list("changepoint.mv.mrc.class"),function(x,p=MDL(x),group=FALSE,display.variable.names=TRUE,show=TRUE)
{
         # nulling out variables used in ggplot to get the package past CRAN checks
         k<-variable<-value<-cpt<-NULL # note, the value of NULL will not really be used - thsi is a hack to satisfy CRAN checks
  
         if(!is.logical(group))
	 {
	   stop("group must be a logical value (TRUE or FALSE)")	
	 }
	 if(!is.numeric(p))
	 {
	    stop("p must be a positive integer")
	 }
	 if(p <= 0)
	 {
	   stop("p must be a positive integer > 0")
	 }
	 if(p > x@pmax)
	 {
	   stop("p must be a positive integer <= pmax")
	 }


	 width<-1.0
	 height<-1.0
	 size<-1.0
	 uv.cpts<-TRUE
	 mr.uv.cpts<-TRUE
	 mr.cpts<-TRUE
	 cpt.colours<-c("yellow","green","red")
	 cpt.alphas<-c(1.0,1.0,1.0)
	 plot.range<-NULL

	 mrc<-x

 n<-nrow(mrc@data)
    
    idx.var<-1:(n)

    df<-mrc@data #  see below


    
variable.names<-names(df)[-1]
# determine the group order

#if(group) {df<-df[,c(1,group.order+1)]}

# add the univariate change points
df.uv.cpts<-Reduce(cbind,Map(function(x) {res<-rep(FALSE,n);res[x]<-TRUE;data.frame(res);},mrc@cpts.uv))
names(df.uv.cpts)<-paste("uv.cpts.",variable.names,sep="")
#if(group) {df.uv.cpts<-df.uv.cpts[,group.order]}    
# add the most recent univariate change points
df.mr.uv.cpts<-Reduce(cbind,Map(function(i) {res<-rep(FALSE,n);res[tail(which(df.uv.cpts[,i]==TRUE),n=2)]<-TRUE;data.frame(res);},1:ncol(df.uv.cpts)))
names(df.mr.uv.cpts)<-paste("mr.uv.cpts.",variable.names,sep="")
#if(group) {df.mr.uv.cpts<-df.mr.uv.cpts[,group.order]}    
                                        # add the most recent changepoints
    
mrcpts<-unlist(Map(function(var) var[length(var)-1],cpts.mv(x,p)))


    
    
df.mr.cpts<-Reduce(cbind,Map(function(i) {res<-rep(FALSE,n);res[i]<-TRUE;data.frame(res);},mrcpts))


    
# df.mr.cpts<-Reduce(cbind,Map(function(i) {res<-rep(1,res@n);res[i]<-2;data.frame(res);},res@mrcpts))
names(df.mr.cpts)<-paste("mr.cpts.",display.variable.names,sep="")
#if(group) {df.mr.cpts<-df.mr.cpts[,group.order]}


# done reversing changepoints
group.order<-NULL    
if(group)
{
   group.order<-order(mrcpts)
}
else
{
   # ggplot is upside down !!
   group.order<-seq(ncol(df)-1,1,-1)     
}
   df<-df[,c(1,group.order+1)]
   df.uv.cpts<-df.uv.cpts[,group.order]
   df.mr.uv.cpts<-df.mr.uv.cpts[,group.order]
   df.mr.cpts<-df.mr.cpts[,group.order]



    
df<-cbind(df,df.uv.cpts,df.mr.uv.cpts,df.mr.cpts)


    
    if(uv.cpts==FALSE) {cpt.alphas[1]=0.0}
    if(mr.uv.cpts==FALSE) {cpt.alphas[2]=0.0}
    if(mr.cpts==FALSE) {cpt.alphas[3]=0.0}
    n.v<-(ncol(df))%/%4
    # scale the data
    df[,2:(n.v+1)]<-apply(df[,2:(n.v+1)],2,function(x) (x-min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE)-min(x,na.rm=TRUE)))
    df.heatmap<-melt(df,id.var="k",measure.vars=names(df)[2:(n.v+1)],value.name="value")
    df.heatmap<-cbind(df.heatmap,
    data.frame("cpt"=unlist(Map(function(a,b,c) if(c) 3 else if(b) 2 else if(a) 1 else 0,
                                melt(df,id.var="k",measure.vars=names(df)[(n.v+2):(2*n.v+1)])[,3],
                                melt(df,id.var="k",measure.vars=names(df)[(2*n.v+2):(3*n.v+1)])[,3],
                                melt(df,id.var="k",measure.vars=names(df)[(3*n.v+2):(4*n.v+1)])[,3]
                                )
                           )
               )
    )
    if(!is.null(plot.range))
    {
      df.heatmap<-df.heatmap[Reduce("|",Map(function(x) df.heatmap$k==x,plot.range)),]
    }
    p <- ggplot(subset(df.heatmap,!is.na(value)),aes(k,variable))
                                        # p <- ggplot(df.heatmap,aes(k,variable))
    if(display.variable.names==FALSE)
        {
            p<-p+theme(axis.text.y=element_blank(),axis.title=element_blank())
        }
    p <- p + geom_tile(aes(fill=value,colour=as.factor(cpt)),width=width,height=height,size=size)
    if(Reduce("|",is.na(df.heatmap$value)))
        {
            p<- p + geom_tile(data=subset(df.heatmap,  is.na(value)), aes(colour = NA),linetype = 0, fill = "white", alpha = 0.5)
        }
    # hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')
    # p <- p + scale_fill_gradientn(colours=hm.palette(100))
    p <- p + scale_colour_manual(values=alpha(c("black",cpt.colours[1],cpt.colours[2],cpt.colours[3]),c(0,cpt.alphas[1],cpt.alphas[2],cpt.alphas[3])))
    p <- p + guides(colour=FALSE)
    if(show==TRUE)
    {
       plot(p)
    }
    invisible(p)


})







#' @name summary
#'
#' @docType methods
#'
#' @rdname summary-method
#'
#' @aliases summary,changepoint.mv.mrc.class-method
#'
# '@examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-spot(mrcexample)
#' summary(res)
#'
#' @export
setMethod("summary",signature=list("changepoint.mv.mrc.class"),function(object)
{
    
    if(length(object@cpts.uv)>0 && length(object@cpts.mv) > 0)
    {
        cat("Univariate and Multivariate Analysis","\n")
    }
    if(length(object@cpts.uv)>0 && length(object@cpts.mv) == 0)
    {
        cat("Univariate Analysis","\n")
    }
    if(length(object@cpts.uv) == 0 && length(object@cpts.mv) > 0)
    {
        cat("Multivariate Analysis","\n")
    }
    #
    #if(length(object@cpts.uv)>0)
    #{
    #    cat("Univariate Changepoints :","\n")
    #    Map(function(varname,cpts) cat(varname," : ",cpts,"\n"),names(object@data)[-1],object@cpts.uv)
    #}
    #if(length(object@cpts.mv)>0)
    #{
    #    cat("Multivariate Changepoints :","\n")
    #    Map(function(varname,cpts) cat(varname," : ",cpts,"\n"),names(object@data)[-1],object@cpts.mv)
    #}
    cat("Cost Function : ",object@cost,"\n")
    #cat("Penalty : ",object@penalty,"\n")
    cat("Alpha Value : ",object@alpha,"\n")
    cat("MDL : ",MDL(object),"\n")
    mat.penalised.costs<-as.matrix(data.frame("p"=1:length(object@penalised.costs),"cost"=object@penalised.costs))
    cat("Penalised Costs : ","\n")
    print(costs(object))
    invisible()
})
















