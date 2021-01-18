#' @rdname predict.evmOpt
#' @export
plot.lp.evmOpt <- function(x, main=NULL,
         pch= 1, ptcol =2 , cex=.75, linecol = 4 ,
         cicol = 1, polycol = 15, plot.=TRUE, ...){
  x <- x$obj
  family <- x$family
  x <- x$link
  
  if(dim(x)[1] == 1){
    stop("Need range of covariate values to plot linear predictors")
  }
  if(length(colnames(x)) <= (length(family$param) + 1) ){
    stop("Please use ci.fit=TRUE in call to predict, to calculate confidence intervals")
  }
  
  makelp <- function(x, family){
      p <- names(family$param)
      res <- vector("list", length=length(p))
      names(res) <- p
      for (i in 1:length(p)){
          res[[i]] <- as.matrix(x[, paste0(p[i], c("", ".lo", ".hi"))])
      }
      res
  }
  Ests <- makelp(data.frame(x), family)
  Names <- names(family$param)
  cn <- colnames(x)
  ParNames <- paste(rep(Names,each=4),c("",".lo",".hi",".se"),sep="")
  which <- ! (cn %in% ParNames)

  X <- x[,which]
  if(is.null(dim(X))){
     X <- matrix(X)
     dimnames(X) <- list(dimnames(x)[[1]],dimnames(x)[[2]][which])
  }
  nPar <- length(Names)
  Output <- list(NULL)
  nPlot <- 0
  for(i in 1:nPar){
    for(j in 1:dim(X)[2]){
      if(length(unique(Ests[[i]][,1])) > 1){
          if(length(unique(X[,j])) > 1){
          ord <- order(X[,j])
          x <- X[ord,j]
          y <- Ests[[i]][ord,]
          nPlot <- nPlot + 1
          Output[[nPlot]] <- list(x=x,y=y,CovName = colnames(X)[i],ParName = Names[i])
          
          if(plot.){
            plot(x, y[,1],type="n",ylab=Names[i],xlab=colnames(X)[j],main=main,ylim=range(y))
              
            if (polycol != 0){
              polygon(c( x,        rev(x)),
                      c(y[,2],rev(y[,3])),
                      col=polycol, border = FALSE) # Close polygon
            } else {
              lines(x, y[,2], col = cicol)
              lines(x, y[,3], col = cicol)
            }

            lines(x, y[,1], col = linecol[ 1 ] )
          }
        }
      }
    }
  }
  invisible(Output)
}

#' @export
plot.lp.evmSim <- function(x, type="median", plot.=TRUE,...){
  if(dim(x$obj$link)[1] == 1){
    stop("Need range of covariate values to plot linear predictors")
  }
  p <- names(x$obj$family$param)
  np <- length(p)
# re-format to same column structure as lp.evmOpt x
  ColIndexMeans <- 1+4*(0:(np-1))
  if(casefold(type) == "median"){
    offset <- 1
  } else if(casefold(type) == "mean") {
    offset <- 0
  } else {
    stop("type must be \"mean\" or \"median\" ")
  }
  which <- c(ColIndexMeans + offset,rep(1:2,np) + rep(ColIndexMeans+1,each=2), (4*np+1): dim(x$obj$link)[2])
  x$obj$link <- x$obj$link[,which]
  colnames(x$obj$link)[1:(3*np)] <-  c(p,paste(rep(p,each=2),rep(c(".lo",".hi"),np),sep=""))

  plot.lp.evmOpt(x,plot.=plot.,...)
}

#' @export
plot.lp.evmBoot <- plot.lp.evmSim

