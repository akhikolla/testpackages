#' @rdname migpd
#' @export
`ggplot.migpd` <-
function(data, mapping=NULL, 
         main=c("Probability plot","Quantile plot","Return level plot","Histogram and density"), 
         xlab=rep(NULL,4), nsim=1000, alpha=.05, ... , environment)
{
    if ( !missing( main ) ){
      if ( length( main ) != 1 & length( main ) != 4 ){
        stop( "main should have length 1 or 4" )
      } else if ( length( main ) == 1 ){ 
        main <- rep( main, 4 ) 
      }
    }
    
    p <- lapply(1:length(data$models), function(i)ggplot(data$model[[i]], main= paste(rep(names(data$model[i]),4),main), xlab=xlab,nsim=nsim,alpha=0.05,plot.=FALSE,..., environment))

    invisible(p)
}


