
### it is suggested to close all graphical devices before plotting DMA results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.grid.dma <- function(x, non.interactive=NULL, ...)
  {

if (is.null(non.interactive)) 
  {
    non.interactive <- FALSE
  }
 
plot1g <- function(x)
  {
    if ( ((ncol(x$RMSE)>=2) && (nrow(x$RMSE)>=2)) && is.na(charmatch("c",paste(rownames(x$RMSE), collapse = ""))) )
      {
        d <- apply(x$RMSE, 2, rev)
        d <- apply(d, 1, rev)

        filled.contour(as.numeric(rownames(d)),as.numeric(colnames(d)),d,col=colorRampPalette(c("blue","green"))(30),
        xlab=expression(paste(alpha)),ylab=expression(paste(lambda)))
         
        rm(d)
      }
  }
 
plot2g <- function(x)
  {
    if ( ((ncol(x$RMSE)>=2) && (nrow(x$RMSE)>=2)) && is.na(charmatch("c",paste(rownames(x$RMSE), collapse = ""))) )
      {
        d <- apply(x$MAE, 2, rev)
        d <- apply(d, 1, rev)

        filled.contour(as.numeric(rownames(d)),as.numeric(colnames(d)),d,col=colorRampPalette(c("blue","green"))(30),
        xlab=expression(paste(alpha)),ylab=expression(paste(lambda)))
         
        rm(d)
      }
  }
 
plot3g <- function(x)
  {
    col <- rich.colors((length(x$models)*length(x$models[[1]])), palette="temperature", rgb=FALSE, plot=FALSE) 
    
    names <- colnames(x[[1]][[1]][[1]]$post.incl)
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x[[1]][[1]][[1]]$post.incl)/7)
      }
    labs <- rownames(x[[1]][[1]][[1]]$post.incl)[inc]
    
    width <-  480
    height <- 300

    for (j in 1:(length(names)))
      {
        p <- 1
        mypath <- file.path(getwd(), paste("grid_pip_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x[[1]][[1]][[1]]$post.incl[,1]), rep(NA,length(index(x[[1]][[1]][[1]]$post.incl[,1]))), lty=2, type="l", col=col[p], ylim=c(0,1), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:length(x$models))
          {
            for (jj in 1:length(x$models[[1]]))
              {
                plot(index(x[[1]][[1]][[1]]$post.incl[,1]),x[[1]][[ii]][[jj]]$post.incl[,j],col=col[p],ylim=c(0,1),
                     axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
                par(new=TRUE)
                p <- p + 1
              }
          }
        axis(1, at=inc, labels=labs)
        dev.off()
      }

     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("grid_pip_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="grid_pip_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
      for(i in 1:(length(names))) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }

plot4g <- function(x)
  {
    col <- rich.colors((length(x$models)*length(x$models[[1]])), palette="temperature", rgb=FALSE, plot=FALSE) 

    names <- colnames(x[[1]][[1]][[1]]$exp.coef.)
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * nrow(x[[1]][[1]][[1]]$exp.coef.)/7)
      }
    labs <- rownames(x[[1]][[1]][[1]]$exp.coef.)[inc]
    
    width <-  480
    height <- 300

    m1 <- vector()
    m2 <- vector()
    
    for (j in 1:(length(names))) 
      {
        m1.t <- x[[1]][[1]][[1]]$exp.coef.[,j]
        for (ii in 1:length(x$models))
          {
            for (jj in 1:length(x$models[[1]]))
              {
                m1.t <- cbind(m1.t, x[[1]][[ii]][[jj]]$exp.coef.[,j])
              }
          }
       m1[j] <- min(m1.t)
       m2[j] <- max(m1.t)
     }

    for (j in 1:(length(names)))  
      {
        p <- 1
        mypath <- file.path(getwd(), paste("grid_coef_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd = TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x[[1]][[1]][[1]]$exp.coef.[,1]), rep(NA,length(index(x[[1]][[1]][[1]]$exp.coef.[,1]))), lty=2, type="l", col=col[p], ylim=c(m1[j],m2[j]), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:length(x$models))
          {
            for (jj in 1:length(x$models[[1]]))
              {
                plot(index(x[[1]][[1]][[1]]$exp.coef.[,1]),x[[1]][[ii]][[jj]]$exp.coef.[,j],col=col[p],ylim=c(m1[j],m2[j]),
                     axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
                par(new=TRUE)
                p <- p + 1
              }
          }
        axis(1, at=inc, labels=labs)
        dev.off()
      }

     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("grid_coef_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="grid_coef_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
      for(i in 1:(length(names))) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }
  
if (x$models[[1]][[1]]$parameters[4]=="DMA")
      {
        if (non.interactive == FALSE) 
          {
            choices <- c("RMSE", "MAE","posterior inclusion probabilities - separate plots (files in working directory)",
                         "expected coefficients - separate plots (files in working directory)" )
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
else
      {
        if (non.interactive == FALSE) 
          {
            choices <- c("RMSE", "MAE",
                         "expected coefficients - separate plots (files in working directory)" )
            pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
            switch(pick, plot1g(x), plot2g(x), plot4g(x))
          }
        else
          {
            plot1g(x)
            plot2g(x)
            plot4g(x)
          }
      }

 
  }
  