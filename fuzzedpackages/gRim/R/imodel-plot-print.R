### #####################################################
###
### Plot iModels
### "D"iscrete variables are "D"ots (=gray)
### "C"ontinuous variabesl are "C"ircles (=transparent)
###
### #####################################################

#' @method iplot iModel
#' @export
iplot.iModel <- function(x,...){
  ig <- ugList(.glist(x), result="igraph")
  V(ig)$label <- V(ig)$name
  V(ig)$size  <- 50
  ig$cex      <-  4
  V(ig)$label.cex <- 1.2

  switch(class(x)[1],
         "dModel"={
           V(ig)$color <- "grey"
         },
         "cModel"={
           V(ig)$color <- "white"
         },
         "mModel"={
           V(ig)$color <- "white"
           disc.idx <- match(datainfo(x)$disc.names, V(ig)$name) #-1
           V(ig)[disc.idx]$color <- "grey"
         })

  ig$layout <- layout.lgl
  plot(ig)
  return(invisible(x))
}

#' @method plot iModel
#' @export
plot.iModel <- function(x,...){
  uG <- ugList(.glist(x))
  switch(class(x)[1],
         "dModel"={
           fillv <- rep("lightgray", length(x$varNames))
           names(fillv) <- x$varNames
         },
         "cModel"={
           fillv <- rep("transparent", length(x$varNames))
           names(fillv) <- x$varNames
         },
         "mModel"={
           dv <- datainfo(x)$disc.names
           cv <- datainfo(x)$cont.names
           fillv  <- c(rep("lightgray", length(dv)), rep("transparent", length(cv)))
           names(fillv) <- c(dv,cv)
         })

  plot(uG, nodeAttrs=list(fillcolor=fillv))
}


#' @method print iModel
#' @export
print.iModel <- function(x, ...){

  cat(sprintf("Model: A %s with %i variables\n", class(x)[1], length(x$varNames)))
  #str(x$varNames)
  ## Model properties
  cat(sprintf(" graphical : %5s  decomposable : %5s\n", x$isGraphical, x$isDecomposable))

  if (x$isFitted){
    dimension <- fitinfo(x)$dimension
    #cat("Fit info: \n")
    cat(sprintf(" -2logL    : %14.2f mdim : %4d aic : %12.2f \n",
                -2*fitinfo(x)$logL,       dimension["mod.dim"], fitinfo(x)$aic))
    cat(sprintf(" ideviance : %14.2f idf  : %4d bic : %12.2f \n",
                fitinfo(x)$ideviance,  dimension["idf"], fitinfo(x)$bic))
    cat(sprintf(" deviance  : %14.2f df   : %4d \n",
                fitinfo(x)$dev,        dimension["df"]))
  }
    
  return(invisible(x))
}

#' @method print dModel
#' @export
print.dModel <- function(x, ...){

  print.iModel(x)
  
  ## If the model is fitted
  ##
  if (x$isFitted){    
    ## Print warnings about sparsity and adjusments of df's
    ##
    if ( (!fitinfo(x)$sparseinfo["df.ok"]) | (!fitinfo(x)$sparseinfo["chi2.ok"])) {
      cat("Notice: Table is sparse\n")
      if (!fitinfo(x)$sparseinfo["chi2.ok"])
        cat(sprintf("  Asymptotic chi2 distribution may be questionable.\n"))
      
      if (!fitinfo(x)$sparseinfo["df.ok"])
        cat(sprintf("  Degrees of freedom can not be trusted.\n"))
      
      if (fitinfo(x)$sparseinfo["sparse.df.ok"] & !fitinfo(x)$sparseinfo["df.ok"]){
        cat(sprintf("  Model dimension adjusted for sparsity : %d\n",
                    fitinfo(x)$dimension["mod.dim.adj"]))
      }
    }
  }
  return(invisible(x))
}



































