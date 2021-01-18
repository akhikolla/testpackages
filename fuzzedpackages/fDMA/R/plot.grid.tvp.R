
### it is suggested to close all graphical devices before plotting TVP results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.grid.tvp <- function(x, non.interactive=NULL, ...)
  {

if (is.null(non.interactive)) 
  {
    non.interactive <- FALSE
  }
   
plot1g <- function(x)
  {
    plot(rownames(x$fq), x$fq[,1], lty=1, type="l", col="blue", ylim=c(min(x$fq[,1]),max(x$fq[,1])), 
         axes=TRUE, xlab=expression(lambda), main="RMSE")
  }
 
plot2g <- function(x)
  {
    plot(rownames(x$fq), x$fq[,2], lty=1, type="l", col="blue", ylim=c(min(x$fq[,2]),max(x$fq[,2])), 
         axes=TRUE, xlab=expression(lambda), main="MAE")
  }
 
plot3g <- function(x)
  {
    col <- rich.colors(nrow(x$fq)+1, palette="temperature", rgb=FALSE, plot=FALSE) 
    
    names <- colnames(x$models[[1]]$thetas)
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x$models[[1]]$y)/7)
      }
    labs <- rownames(x$models[[1]]$y)[inc]
    
    width <-  480
    height <- 300

    m1 <- vector()
    m2 <- vector()
    
    for (j in 1:(length(names))) 
      {
        m1.t <- x$models[[1]]$thetas[,j]
        for (ii in 1:length(x$models))
          {
            m1.t <- cbind(m1.t, x$models[[ii]]$thetas[,j])
          }
       m1[j] <- min(m1.t)
       m2[j] <- max(m1.t)
     }

    for (j in 1:(length(names)))  
      {
        mypath <- file.path(getwd(), paste("grid_tvp_coef_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x$models[[1]]$y), rep(NA,length(index(x$models[[1]]$y))), lty=2, type="l", col=col[1], ylim=c(m1[j],m2[j]), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:length(x$models))
          {
            plot(index(x$models[[1]]$y),x$models[[ii]]$thetas[,j],col=col[ii+1],ylim=c(m1[j],m2[j]),
                 axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
            par(new=TRUE)
          }
        axis(1, at=inc, labels=labs)
        dev.off()
      }

     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("grid_tvp_coef_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="grid_tvp_coef_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
      for(i in 1:(length(names))) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }

        if (non.interactive == FALSE) 
          {
            choices <- c("RMSE", "MAE",
                         "coefficients - separate plots (files in working directory)" )
            pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
            switch(pick, plot1g(x), plot2g(x), plot3g(x))
          }
        else
          {
            plot1g(x)
            plot2g(x)
            plot3g(x)
          }
 
  }
  