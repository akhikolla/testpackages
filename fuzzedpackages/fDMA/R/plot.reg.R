
### it is suggested to close all graphical devices before plotting "reg" class results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.reg <- function(x, non.interactive=NULL, ...)
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
    col <- rich.colors(ncol(x$coeff.), palette="temperature", rgb=FALSE, plot=FALSE)

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]

    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1)) 
    
    m1 <- min(x$coeff.,na.rm=TRUE)
    m2 <- max(x$coeff.,na.rm=TRUE)
    
    for (i in 1:ncol(x$coeff.))
      {
        if (i==ncol(x$coeff.))
          { 
            axes <- TRUE
            main <- "regression coefficients"
          }
        else
          {
            axes <- FALSE
            main <- ""
          }
        plot(x$coeff.[,i], type="l", col=col[i], ylim=c(m1,m2), axes=axes, xaxt='n', xlab="", ylab="", main=main)
        if (! i==ncol(x$coeff.)) { par(new=TRUE) }
      }
    axis(1, at=inc, labels=labs)
    legend('bottom', inset=c(0,-0.5), colnames(x$coeff.), lty=rep(1,i), col=col[1:i], ncol=6, cex=0.6) 
  }

plot4 <- function(x)
  {
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$coeff.)/7)
      }
    labs <- rownames(x$y)[inc]
    names <- colnames(x$coeff.)
    
    width <-  480
    height <- 300

    for (j in 1:(ncol(x$coeff.)))
      {
        mypath <- file.path(getwd(), paste("reg_coeff_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        m1 <- min(x$coeff.[,j],na.rm=TRUE)
        m2 <- max(x$coeff.[,j],na.rm=TRUE)
        plot(x$coeff.[,j],col="blue",ylim=c(m1,m2),axes=TRUE, xaxt='n', xlab='',ylab='',type="l",main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
       }

     img <- list()
     for (i in 1:j)
      {
        mypath <- file.path(getwd(), paste("reg_coeff_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="reg_coeff_all.png", width = 2 * width, height = height * ceiling(j/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling(j/2)), ncol=2, byrow=TRUE))
      for(i in 1:j) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }

plot5 <- function(x)
  {
    col <- rich.colors(ncol(x$coeff.), palette="temperature", rgb=FALSE, plot=FALSE)

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]

    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1)) 
    
    for (i in 1:ncol(x$coeff.))
      {
        if (i==ncol(x$coeff.))
          { 
            axes <- TRUE
            main <- "p-values for t-test"
          }
        else
          {
            axes <- FALSE
            main <- ""
          }
        plot(x$p.val.[,i], type="l", col=col[i], ylim=c(0,1), axes=axes, xaxt='n', xlab="", ylab="", main=main)
        if (! i==ncol(x$coeff.)) { par(new=TRUE) }
      }
    axis(1, at=inc, labels=labs)
    legend('bottom', inset=c(0,-0.5), colnames(x$coeff.), lty=rep(1,i), col=col[1:i], ncol=6, cex=0.6) 
  }

plot6 <- function(x)
  {
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$coeff.)/7)
      }
    labs <- rownames(x$y)[inc]
    names <- colnames(x$coeff.)
    
    width <-  480
    height <- 300

    for (j in 1:(ncol(x$coeff.)))
      {
        mypath <- file.path(getwd(), paste("reg_p_val_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(x$p.val.[,j],col="blue",ylim=c(0,1),axes=TRUE, xaxt='n', xlab='',ylab='',type="l",main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
       }

     img <- list()
     for (i in 1:j)
      {
        mypath <- file.path(getwd(), paste("reg_p_val_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="reg_p_val_all.png", width = 2 * width, height = height * ceiling(j/2))
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
    choices <- c("actual and predicted", "residuals", 
                 "coefficients - one plot", "coefficients - separate plots (files in working directory)",
                 "p-values - one plot", "p-values - separate plots (files in working directory)")
    pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
    switch(pick, plot1(x), plot2(x), plot3(x), plot4(x), plot5(x), plot6(x))
  }
else
  {
    plot1(x)
    plot2(x)
    plot3(x)
    plot4(x)
    plot5(x)
    plot6(x)  
  }

  }
  
  