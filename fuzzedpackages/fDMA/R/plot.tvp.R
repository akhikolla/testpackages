
### it is suggested to close all graphical devices before plotting "tvp" class results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.tvp <- function(x, non.interactive=NULL, ...)
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
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
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
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
    resid <- as.numeric(x$y - x$y.hat)
    par(xpd=FALSE)
    plot(resid, type="l", col="blue", ylim=c(min(resid),max(resid)), axes=TRUE, xaxt='n', xlab="", ylab="", main="residuals")
    axis(1, at=inc, labels=labs)
  }

plot3 <- function(x)
  {
    col <- rich.colors(ncol(x$thetas), palette="temperature", rgb=FALSE, plot=FALSE)

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(as.vector(x$y))/7)
      }
    labs <- rownames(x$y)[inc]

    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1)) 
    
    m1 <- min(x$thetas,na.rm=TRUE)
    m2 <- max(x$thetas,na.rm=TRUE)
    
    for (i in 1:ncol(x$thetas))
      {
        if (i==ncol(x$thetas))
          { 
            axes <- TRUE
            main <- "regression coefficients"
          }
        else
          {
            axes <- FALSE
            main <- ""
          }
        plot(x$thetas[,i], type="l", col=col[i], ylim=c(m1,m2), axes=axes, xaxt='n', xlab="", ylab="", main=main)
        if (! i==ncol(x$thetas)) { par(new=TRUE) }
      }
    axis(1, at=inc, labels=labs)
    legend('bottom', inset=c(0,-0.5), colnames(x$thetas), lty=rep(1,i), col=col[1:i], ncol=6, cex=0.6) 
  }

plot4 <- function(x)
  {
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$thetas)/7)
      }
    labs <- rownames(x$y)[inc]
    names <- colnames(x$thetas)
    
    width <-  480
    height <- 300

    for (j in 1:(ncol(x$thetas)))
      {
        mypath <- file.path(getwd(), paste("tvp_coeff_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        m1 <- min(x$thetas[,j],na.rm=TRUE)
        m2 <- max(x$thetas[,j],na.rm=TRUE)
        plot(x$thetas[,j],col="blue",ylim=c(m1,m2),axes=TRUE, xaxt='n', xlab='',ylab='',type="l",main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
       }

     img <- list()
     for (i in 1:j)
      {
        mypath <- file.path(getwd(), paste("tvp_coeff_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="tvp_coeff_all.png", width = 2 * width, height = height * ceiling(j/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling(j/2)), ncol=2, byrow=TRUE))
      for(i in 1:j) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }

if (non.interactive == FALSE) 
  {
    choices <- c("actual and predicted", "residuals", "coefficients - one plot", "coefficients - separate plots (files in working directory)")
    pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
    switch(pick, plot1(x), plot2(x), plot3(x), plot4(x))
  }
else
  {
    plot1(x)
    plot2(x)
    plot3(x)
    plot4(x)  
  }
  

  }
  
  