
### it is suggested to close all graphical devices before plotting altf3() results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.altf3 <- function(x, non.interactive=NULL, ...)
  {
 
plot1g <- function(x)
  {
  
    if (is.null(non.interactive)) 
      {
        non.interactive <- FALSE
      }
    
    names <- colnames(x$coeff.[[1]])
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$coeff.[[1]])/7)
      }
    labs <- rownames(x$y)[inc]
    
    width <-  480
    height <- 300
    

    for (j in 1:(length(names)))
      {
        m1 <- min(x$coeff.[[1]][,j],na.rm=TRUE)
        m2 <- max(x$coeff.[[1]][,j],na.rm=TRUE)

        mypath <- file.path(getwd(), paste("altf3_coeff_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x$y), x$coeff.[[1]][,j], lty=1, type="l", col="blue", ylim=c(m1,m2), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
      }

     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("altf3_coeff_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="altf3_coeff_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
      for(i in 1:(length(names))) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }


plot2g <- function(x)
  {
    
    names <- colnames(x$p.val.[[1]])
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$p.val.[[1]])/7)
      }
    labs <- rownames(x$y)[inc]
    
    width <-  480
    height <- 300
    

    for (j in 1:(length(names)))
      {
        m1 <- 0
        m2 <- 1

        mypath <- file.path(getwd(), paste("altf3_p_val_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x$y), x$p.val.[[1]][,j], lty=1, type="l", col="blue", ylim=c(m1,m2), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        axis(1, at=inc, labels=labs)
        dev.off()
      }

     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("altf3_p_val_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="altf3_p_val_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
      for(i in 1:(length(names))) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }


plot3g <- function(x)
  {
    col <- rich.colors(ncol(x$weights[[1]]), palette="temperature", rgb=FALSE, plot=FALSE) 
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$weights[[1]])/7)
      }
    labs <- rownames(x$y)[inc]
    
    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1))
  
    for (i in 1:(ncol(x$weights[[1]])-1))
      {
        plot(x$weights[[1]][,i], type="l", col=col[i], ylim=c(0,1), axes=FALSE, xaxt='n', xlab="", ylab="", main="")
        par(new=TRUE)
      }
    plot(x$weights[[1]][,i+1], type="l", col=col[i+1], ylim=c(0,1), axes=TRUE, xaxt='n', xlab="", ylab="", main="models' weights")
    axis(1, at=inc, labels=labs)
    legend('bottom', inset=c(0,-0.45), colnames(x$weights[[1]]), lty=rep(1,(i+1)), col=col[1:(i+1)], ncol=5, cex=0.6) 

  }


plot4g <- function(x)
  {

    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$weights[[1]])/7)
      }
    labs <- rownames(x$y)[inc]
    
    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
    
    windows <- as.numeric(colnames(x$weights[[1]]))
    
    plot(as.vector(x$exp.win.[[1]]), type="l", col="blue", ylim=c(min(windows),max(windows)), axes=TRUE, xaxt='n', xlab="", ylab="", main="exp. window size")
    axis(1, at=inc, labels=labs)

  }

        if (non.interactive == FALSE) 
          {
            choices <- c("expected coefficients - separate plots (files in working directory)",
                         "p-values for t-tests - separate plots (files in working directory)",
                         "models' weights",
                         "expected window size")
            pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
            switch(pick, plot1g(x), plot2g(x), plot3g(x), plot4g(x))
          }
        else
          {
            plot1g(x)
            plot2g(x)
            plot3g(x)
            plot4g(x)
          }
 
  }
  