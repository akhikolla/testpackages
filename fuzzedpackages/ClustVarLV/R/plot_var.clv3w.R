#' Scores plot from a Candecomp Parafac analysis. The group membership of the variables is superimposed.
#'
#' @param resclv3w the data matrix
#' @param K  the number of groups in the partition (already defined if CLV3W_kmeans i used)
#' @param axeh  component number for the horizontal axis
#' @param axev  component number for the vertical axis
#' @param labels boolean to add variable' labels (label=TRUE) on the plot or not (label=FALSE). By default label=TRUE
#' @param cex.lab magnification to be used for labels (1 by default)
#' @param v_colors default NULL. If missing colors are given, by default
#' @param v_symbol symbols are given instead of colors for the identification of the groups/ =FALSE: no symbol (by default).
#' @param beside plot per cluster of variables, side-by-side\cr
#'  =FALSE : an unique plot with all the variables with the identification of their group membership (by default).
#' @param mode3 projection of the mode 3 elements onto the scores plot\cr
#'  =FALSE : mode 3 elements are not represented (by default).
#'
#' @export
#'

plot_var.clv3w <- function(resclv3w,K=NULL,axeh = 1, axev = 2, labels=FALSE,
                  cex.lab = 1, v_colors = NULL, v_symbol = FALSE, beside = FALSE,mode3=FALSE) {

  if (!inherits(resclv3w, "clv3w") )   stop("non convenient object")
  appel      <- as.list(resclv3w$call)
  X          <- resclv3w$param$X
  X <- block.scale(X, xcenter = TRUE, xscale = appel$mode.scale)
  n <- dim(X)[[1]]
  p <- dim(X)[[2]]
  q <- dim(X)[[3]]
  if (axeh > axev) {
    temp=axeh
    axeh <- axev
    axev <- temp
  }
    if (is.null(v_colors)) {
      v_colors <- c("blue", "red", "green", "black", "purple",
                    "orange", "yellow", "tomato", "pink", "gold", "cyan",
                    "turquoise", "violet", "green", "khaki", "maroon",
                    "chocolate", "burlywood")
    }
    if (v_symbol) {
      v_colors <- rep("black", 20)
    }

  # candecomp parafac
    cp2 <- CP2_MS(X,ncp=axev)
    I.tot <- cp2$afit[1]+cp2$aloss[1]
    cp.in1 <- round(cp2$afit[axeh]/I.tot*100,2)
    cp.in2 <- round(cp2$afit[axev]/I.tot*100,2)
    scores <- cbind(cp2$u[,axeh]*sqrt(cp2$afit[axeh]),cp2$u[,axev]*sqrt(cp2$afit[axev]))
    vp<-(scores[,1]^2+scores[,2]^2)

    coordvar <- cbind(cp2$v[,axeh],cp2$v[,axev])
    coordind <- scores
    par(pty="s")

    lp=max(vp)
    vv<-(coordvar[,1]^2+coordvar[,2]^2)
    lv=max(vv)
    f=sqrt(lp/lv)
    par(pty="s")
    posi <- rep(1,length=n)
    plot(c(coordvar[,1]*f,coordind[,1]),c(coordvar[,2]*f,coordind[,2]),type="n",
         xlab=paste("Dim ",1,"(",round(cp.in1,2),"%)"),
         ylab=paste("Dim ",2,"(",round(cp.in2,2),"%)"),
         main="Candecomp Parafac Scores plot")
    text(scores[,1],scores[,2],labels=dimnames(X)[[1]],cex=cex.lab)
    abline(h=0,v=0)
    par(pty="m")
    if(is.null(eval.parent(appel$K))) {
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      if ((K >0) & (K<resclv3w$param$gmax+1)) {
        clusters <- resclv3w[[K]]$clusters[2,]
        comp     <- resclv3w[[K]]$comp
        loading  <- resclv3w[[K]]$loading
        weight   <- round(resclv3w[[K]]$weight,3)
        criterion<- resclv3w[[K]]$criterion
        pk <- table(clusters)
      }
      else stop("invalid number of clusters")
    } else {
      clusters<-resclv3w$clusters[2,]
      K<-eval.parent(appel$K)
      comp     <- resclv3w$comp
      loading  <- resclv3w$loading
      weight   <- round(resclv3w$weight,3)
      criterion<- resclv3w$criterion
      pk <- table(clusters)
    }
    colpart <- NULL
    symbpart <- NULL
    symbpart <- 20
    for (j in 1:p) {
      colpart[j] <- v_colors[clusters[j]]
      if (v_symbol)
        symbpart[j] = clusters[j]
    }

    color.var <-colpart
    if (beside==FALSE) {
      dev.new()
      par(pty="s")

      #plot slice 2
      plot(c(coordvar[,1]*f,coordind[,1]),c(coordvar[,2]*f,coordind[,2]),type="n",
           xlab=paste("Dim ",1,"(",round(cp.in1,2),"%)"),
           ylab=paste("Dim ",2,"(",round(cp.in2,2),"%)"),
           main="Candecomp Parafac Scores plot")


      legend("topright",col=v_colors,legend=paste("Cluster",1:K),lty=1, cex=0.8,box.lty=0)

      arrows(0,0,coordvar[,1]*f,coordvar[,2]*f,length=0.1,angle=10,lwd=0.5,col=color.var)
      if (labels==TRUE)
        text(coordvar[,1]*f,coordvar[,2]*f,labels=dimnames(X)[[2]],col=color.var)
      else
        points(coordvar[,1]*f,coordvar[,2]*f,pch=symbpart)
      posi=rep(1,n)
      posi[which(coordind[,1]>max(c(coordvar[,1]*f,coordind[,1]))*0.8)]=2
      posi[which(coordind[,1]<min(c(coordvar[,1]*f,coordind[,1]))*0.8)]=4
      posi[which(coordind[,2]>max(c(coordvar[,2]*f,coordind[,2]))*0.8)]=1
      posi[which(coordind[,2]<min(c(coordvar[,2]*f,coordind[,2]))*0.8)]=3
      text(coordind[,1],coordind[,2],labels=dimnames(X)[[1]],pos=posi,cex=cex.lab)
      abline(h=0,v=0)
    }
    else {
      dev.new()
      nr <- round(sqrt(K),0)
      nc <- K %/% nr
      if (K %% nr )
        nc <- nc+1
      par(mfrow = c(nr, nc))
      case = levels(as.factor(colpart))
      for (k in 1:K) {
         who = which(clusters == k)
         #plot slice 2
         plot(c(coordvar[,1]*f,coordind[,1]),c(coordvar[,2]*f,coordind[,2]),type="n",
              xlab=paste("Dim ",1,"(",round(cp.in1,2),"%)"),
              ylab=paste("Dim ",2,"(",round(cp.in2,2),"%)"),
              main="Candecomp Parafac Scores plot")
         legend("topright",col=v_colors[k],legend=paste("Cluster",k),lty=1, cex=0.8,box.lty=0)
         arrows(0,0,coordvar[who,1]*f,coordvar[who,2]*f,length=0.1,angle=10,lwd=0.5,col=color.var[who])
         if (labels==TRUE)
           text(coordvar[who,1]*f,coordvar[who,2]*f,labels=dimnames(X)[[2]][who],col=color.var[who])
         else
           points(coordvar[,1]*f,coordvar[,2]*f,pch=symbpart[who])
         posi=rep(1,n)
         posi[which(coordind[,1]>max(c(coordvar[,1]*f,coordind[,1]))*0.8)]=2
         posi[which(coordind[,1]<min(c(coordvar[,1]*f,coordind[,1]))*0.8)]=4
         posi[which(coordind[,2]>max(c(coordvar[,2]*f,coordind[,2]))*0.8)]=1
         posi[which(coordind[,2]<min(c(coordvar[,2]*f,coordind[,2]))*0.8)]=3
         text(coordind[,1],coordind[,2],labels=dimnames(X)[[1]],pos=posi,cex=cex.lab)
      }
    }

    #plot slice 3
    if (mode3==TRUE) {
      coordvar <- cbind(cp2$w[,axeh],cp2$w[,axev])
      coordind <- scores
      par(pty="s")

      vp<-(coordind[,1]^2+coordind[,2]^2)
      lp=max(vp)
      vv<-(coordvar[,1]^2+coordvar[,2]^2)
      lv=max(vv)
      f=sqrt(lp/lv)
      dev.new()
      plot(c(coordvar[,1]*f,coordind[,1]),c(coordvar[,2]*f,coordind[,2]),type="n",
           xlab=paste("Dim ",1,"(",round(cp.in1,2),"%)"),
           ylab=paste("Dim ",2,"(",round(cp.in2,2),"%)"),
           main="Candecomp Parafac Scores plot")
      palet <- rainbow(q)
      legend("topright",col=palet,legend=dimnames(X)[[3]],lty=1, cex=0.8,box.lty=0)

      arrows(0,0,coordvar[,1]*f,coordvar[,2]*f,length=0.1,angle=10,lwd=0.5,col=palet)
      if (labels==TRUE)
        text(coordvar[,1]*f,coordvar[,2]*f,labels=dimnames(X)[[3]],col=palet)
      posi=rep(1,n)
      posi[which(coordind[,1]>max(c(coordvar[,1]*f,coordind[,1]))*0.8)]=2
      posi[which(coordind[,1]<min(c(coordvar[,1]*f,coordind[,1]))*0.8)]=4
      posi[which(coordind[,2]>max(c(coordvar[,2]*f,coordind[,2]))*0.8)]=1
      posi[which(coordind[,2]<min(c(coordvar[,2]*f,coordind[,2]))*0.8)]=3
      text(coordind[,1],coordind[,2],labels=dimnames(X)[[1]],pos=posi,cex=cex.lab)
    }
  par(pty="m")
  par(mfrow = c(1, 1))
}

