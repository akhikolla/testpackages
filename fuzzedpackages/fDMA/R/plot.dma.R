
### it is suggested to close all graphical devices before plotting DMA results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.dma <- function(x, non.interactive=NULL, ...)
  {

if (is.null(non.interactive)) 
  {
    non.interactive <- FALSE
  }
       
plot1 <- function(x)
  {
  
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]
    par(xpd=TRUE)
    plot(x$y.hat, type="l", col="blue", ylim=c(min(x$y,x$y.hat),max(x$y,x$y.hat)), axes=FALSE, xaxt='n', xlab="", ylab="", main="")
    par(new=TRUE)
    plot(x$y, type="l", col="red", ylim=c(min(x$y,x$y.hat),max(x$y,x$y.hat)), xaxt='n', xlab="", ylab="", main="actual and predicted")
    axis(1, at=inc, labels=labs)
    legend('bottom', inset=c(0,-0.45), c("actual","predicted"), lty=c(1,1), col=c("red","blue"), ncol=2) 
  }

plot2 <- function(x)
  {
  
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]
    resid <- as.numeric(x$y - x$y.hat)
    par(xpd=FALSE)
    plot(resid, type="l", col="blue", ylim=c(min(resid),max(resid)), axes=TRUE, xaxt='n', xlab="", ylab="", main="residuals")
    axis(1, at=inc, labels=labs)
  }

plot3 <- function(x)
  {
  
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]
    par(xpd=FALSE)
    plot(x$exp.var, type="l", col="blue", ylim=c(0,ncol(x$models)), axes=TRUE, xaxt='n', xlab="", ylab="", main="exp. var")
    axis(1, at=inc, labels=labs)
  }

plot4 <- function(x)
  {
    col <- rich.colors(ncol(x$post.incl), palette="temperature", rgb=FALSE, plot=FALSE)

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]

    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1)) 
    
    if (! x$parameters[1,4] == "DMA") { temp <- replace(x$post.incl,which(x$post.incl==0),NA) }
    
    for (i in 1:(ncol(x$post.incl)-1))
      {
        if (x$parameters[1,4] == "DMA")
          {
            plot(x$post.incl[,i], type="l", col=col[i], ylim=c(0,1), axes=FALSE, xaxt='n', xlab="", ylab="", main="")
          }
        if (x$parameters[1,4] == "DMS" || x$parameters[1,4] == "MED")
          {
            plot(i/ncol(temp)*temp[,i], type="p", pch=15, col=col[i], ylim=c(0,1), axes=FALSE, xaxt='n', xlab="", ylab="", main="")
          }
        par(new=TRUE)
      }
    if (x$parameters[1,4] == "DMA")
      {
        plot(x$post.incl[,i+1], type="l", col=col[i+1], ylim=c(0,1), axes=TRUE, xaxt='n', xlab="", ylab="", main="posterior inclusion probabilities")
      }
    if (x$parameters[1,4] == "DMS" || x$parameters[1,4] == "MED")
      {
        plot(temp[,i+1], type="p", pch=15, col=col[i+1], ylim=c(0,1), axes=TRUE, xaxt='n', yaxt='n', xlab="", ylab="", main="inclusion of variables")
      }
    axis(1, at=inc, labels=labs)
    if (x$parameters[1,4] == "DMA") 
      { 
        legend('bottom', inset=c(0,-0.5), colnames(x$post.incl), lty=rep(1,(i+1)), col=col[1:(i+1)], ncol=6, cex=0.6) 
      }
    if (x$parameters[1,4] == "DMS" || x$parameters[1,4] == "MED")
      {
        legend('bottom', inset=c(0,-0.5), colnames(x$post.incl), pch=rep(c(15),(i+1)), col=col[1:(i+1)], ncol=6, cex=0.6) 
      }
  }

plot5 <- function(x)
  {
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]
    names <- colnames(x$models)
    
    width <-  480
    height <- 300

    for (j in 1:(ncol(x$post.incl)))
      {
        mypath <- file.path(getwd(), paste("pip_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(x$post.incl[,j],col="blue",ylim=c(0,1),axes=TRUE, xaxt='n', xlab='',ylab='',type="l",main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
       }

     img <- list()
     for (i in 1:j)
      {
        mypath <- file.path(getwd(), paste("pip_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="pip_all.png", width = 2 * width, height = height * ceiling(j/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling(j/2)), ncol=2, byrow=TRUE))
      for(i in 1:j) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }

plot6 <- function(x)
  {
    col <- rich.colors(ncol(x$exp.coef.)+1, palette="temperature", rgb=FALSE, plot=FALSE)

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$exp.coef.)/7)
      }
    labs <- rownames(x$exp.coef.)[inc]

    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1)) 
    
    for (i in 1:(ncol(x$exp.coef.)))
      {
        plot(x$exp.coef.[,i], type="l", col=col[i], ylim=c(min(x$exp.coef.),max(x$exp.coef.)), axes=FALSE, xaxt='n', xlab="", ylab="", main="")
        par(new=TRUE)
      }
    plot(index(x$exp.coef.[,i]), rep(0,length(x$exp.coef.[,i])), lty=2, type="l", col=col[i+1], ylim=c(min(x$exp.coef.),max(x$exp.coef.)), 
         axes=TRUE, xaxt='n', xlab="", ylab="", main="coefficients")
    axis(1, at=inc, labels=labs)
    legend('bottom', inset=c(0,-0.5), colnames(x$exp.coef.), lty=rep(1,i), col=col[1:i], ncol=6, cex=0.6) 
  }

plot7 <- function(x)
  {
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$exp.coef.)/7)
      }
    labs <- rownames(x$exp.coef.)[inc]
    names <- colnames(x$models)
    
    width <-  480
    height <- 300

    for (j in 1:(ncol(x$exp.coef.)))
      {
        mypath <- file.path(getwd(), paste("coef_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x$exp.coef.[,j]),rep(0,length(x$exp.coef.[,j])),lty=2,col="black",ylim=c(min(x$exp.coef.[,j]),max(x$exp.coef.[,j])),
             axes=FALSE, xaxt='n', xlab='',ylab='',type="l",main='')
        par(new=TRUE)
        plot(x$exp.coef.[,j],col="blue",ylim=c(min(x$exp.coef.[,j]),max(x$exp.coef.[,j])),axes=TRUE, xaxt='n', xlab='',ylab='',type="l",main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
       }

     img <- list()
     for (i in 1:j)
      {
        mypath <- file.path(getwd(), paste("coef_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="coef_all.png", width = 2 * width, height = height * ceiling(j/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling(j/2)), ncol=2, byrow=TRUE))
      for(i in 1:j) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }
   
plot8 <- function(x)
  {
    col <- rich.colors(ncol(x$post.mod), palette="temperature", rgb=FALSE, plot=FALSE)

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]

    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1)) 
    
    for (i in 1:(ncol(x$post.mod)-1))
      {
        plot(x$post.mod[,i], type="l", col=col[i], ylim=c(0,1), axes=FALSE, xaxt='n', xlab="", ylab="", main="")
        par(new=TRUE)
      }
    plot(x$post.mod[,i+1], type="l", col=col[i+1], ylim=c(0,1), axes=TRUE, xaxt='n', xlab="", ylab="", main="posterior model probabilities")
    axis(1, at=inc, labels=labs)
  }    
     
plot9 <- function(x)
  {
  
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]
    par(xpd=FALSE)
    plot(x$DOW.n.mods.t, type="l", col="blue", ylim=c(0,max(x$DOW.n.mods.t)), axes=TRUE, xaxt='n', xlab="", ylab="", main="number of models used in DMA")
    axis(1, at=inc, labels=labs)
  }
  
plot10 <- function(x)
  {
  
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$post.incl)/7)
      }
    labs <- rownames(x$post.incl)[inc]
    par(xpd=FALSE)
    plot(x$exp.lambda, type="l", col="blue", ylim=c(min(x$exp.lambda),max(x$exp.lambda)), axes=TRUE, xaxt='n', xlab="", ylab="", main="exp. lambda")
    axis(1, at=inc, labels=labs)
  }


    if (x$parameters[1,4] == "DMA")
      {
        if (anyNA(x$DOW.n.mods.t))
          {
            if (non.interactive == FALSE) 
              {
                choices <- c("actual and predicted", "residuals","exp. var", "posterior inclusion probabilities - one plot", 
                             "posterior inclusion probabilities - separate plots (files in working directory)",
                             "expected coefficients - one plot", "expected coefficients - separate plots (files in working directory)", 
                             "exp. lambda", "posterior model probabilities")
                pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
                switch(pick, plot1(x), plot2(x), plot3(x), plot4(x), plot5(x), plot6(x), plot7(x), plot10(x), plot8(x))
              }
            else
              {
                plot1(x)
                plot2(x)
                plot3(x)
                plot4(x)
                plot5(x)
                plot6(x)
                plot7(x)
                plot10(x)
                plot8(x)
              }
          }
        else
          {
            if (non.interactive == FALSE) 
              {
                choices <- c("actual vs. predicted", "residuals","exp. var", "posterior inclusion probabilities - one plot", 
                             "posterior inclusion probabilities - separate plots (files in working directory)",
                             "expected coefficients - one plot", "expected coefficients - separate plots (files in working directory)", 
                             "exp. lambda", "number of models used in DMA estimation")
                pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
                switch(pick, plot1(x), plot2(x), plot3(x), plot4(x), plot5(x), plot6(x), plot7(x), plot10(x), plot9(x))
              }
            else
              {
                plot1(x)
                plot2(x)
                plot3(x)
                plot4(x)
                plot5(x)
                plot6(x)
                plot7(x)
                plot10(x)
                plot9(x)
              }
          }
      }
    
    if (x$parameters[1,4] == "DMS" || x$parameters[1,4] == "MED")
      {
         if (non.interactive == FALSE) 
           {
             choices <- c("actual vs. predicted", "residuals","exp. var", "posterior inclusion probabilities", 
                          "expected coefficients - one plot", "expected coefficients - separate plots (files in working directory)",
                          "exp. lambda")
             pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
             switch(pick, plot1(x), plot2(x), plot3(x), plot4(x), plot6(x), plot7(x), plot10(x))
           }
         else
          {
            plot1(x)
            plot2(x)
            plot3(x)
            plot4(x)
            plot6(x)
            plot7(x)
            plot10(x)
          }
      }

  }
  
  