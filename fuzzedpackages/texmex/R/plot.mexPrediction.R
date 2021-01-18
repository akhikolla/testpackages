#' @rdname mex
#' @export
`plot.predict.mex` <-
    function( x, pch=c( 1, 3, 20 ), col=c( 2, 8, 3), cex=c( 1, 1, 1 ), ask = TRUE, ... ){

    d <- dim( x$data$simulated )[[ 2 ]] -1
    if ( prod( par( "mfrow" ) ) < d ){
        if ( ask ) {
            op <- par(ask = TRUE)
            on.exit(par(op))
        }
    }

  xdat <- x$data$real[, 1 ]
  upts <- seq(from =0.001,to=1-0.0001,len=100)
  xpts <- revTransform(upts,data=x$data$real[, 1 ], qu = mean(x$data$real[,1] < x$mth[1]), th=x$mth[1],sigma = x$gpd.coef[3,1], xi = x$gpd.coef[4,1])

	for( i in 2:( dim( x$data$real )[[ 2 ]] ) ){
		ydat <- x$data$real[, i ]
		xlimits <- range( xdat , x$data$simulated[ , 1 ] )
		ylimits <- range( ydat , x$data$simulated[ , i ] )

		plot( xdat , ydat , xlim=xlimits , ylim=ylimits,
			  xlab = names( x$data$simulated )[ 1 ],
			  ylab = names( x$data$simulated )[ i ],
			  type = "n",...
			 )
		points( x$data$simulated[ x$data$CondLargest, 1 ], x$data$simulated[ x$data$CondLargest, i ], col=col[ 3 ], pch=pch[ 3 ], cex=cex[ 3 ] )
		points( x$data$simulated[!x$data$CondLargest, 1 ], x$data$simulated[!x$data$CondLargest, i ], col=col[ 2 ], pch=pch[ 2 ], cex=cex[ 2 ] )
		points( xdat, ydat , pch=pch[ 1 ], col=col[ 1 ], cex= cex[ 1 ] )
		abline( v = x$data$pth, lty=2, col=3 )
    ypts <- revTransform(upts,data=x$data$real[, i ], qu = mean(x$data$real[,i] < x$mth[i]), th=x$mth[i],sigma = x$gpd.coef[3,i], xi = x$gpd.coef[4,i])
    lines(xpts,ypts,col=3)
	}

	invisible()
}

#' @rdname mex
#' @export
`ggplot.predict.mex` <-
    function(data=NULL, mapping, xlab, ylab,  main,
             ptcol=c("grey","dark blue","orange"), col="dark blue", fill="orange", shape=16:18, size=rep(1,3), plot.=TRUE,..., environment){
        
        xdat <- data$data$real[, 1 ]
        upts <- seq(from =0.001,to=1-0.0001,len=100)
        xpts <- revTransform(upts,data=xdat, qu = mean(xdat < data$mth[1]), th=data$mth[1],
                             sigma = data$gpd.coef[3,1], xi = data$gpd.coef[4,1])
        
        plotfn <- function(i){
            dat <- data.frame(x=xdat,y=data$data$real[, i])
            diag <- data.frame(x=xpts, y=revTransform(upts,data=data$data$real[, i ], qu = mean(data$data$real[,i] < data$mth[i]), th=data$mth[i],sigma = data$gpd.coef[3,i], xi = data$gpd.coef[4,i]))
            cl <- data.frame(x=data$data$simulated[ data$data$CondLargest, 1 ], 
                             y=data$data$simulated[ data$data$CondLargest, i ])
            NOTcl <- data.frame(x=data$data$simulated[!data$data$CondLargest, 1 ], 
                                y=data$data$simulated[!data$data$CondLargest, i ])
    
            ggplot(dat,aes(x,y)) +
                labs(x=names( data$data$simulated )[ 1 ],y=names( data$data$simulated )[ i ]) +
                geom_point(data=cl,shape=shape[2],size=size[2],colour=ptcol[2],alpha=0.5) + 
                geom_point(data=NOTcl,shape=shape[3],size=size[3],colour=ptcol[3],alpha=0.5) + 
                geom_point(shape=shape[1],size=size[1],colour=ptcol[1],alpha=0.5) +
                geom_vline(xintercept=data$data$pth, lty=2, colour=col) +
                geom_line(data=diag,colour=fill)
        }
        
        p <- lapply(2:( dim( data$data$real )[[ 2 ]] ), plotfn)
        rowCol <- grDevices::n2mfrow(length(p))
        
        if (plot.) suppressWarnings(do.call("grid.arrange", c(p, list(nrow=rowCol[1], ncol=rowCol[2]))))
        invisible(p)
    }

