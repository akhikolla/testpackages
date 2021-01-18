
### it is suggested to close all graphical devices before plotting altf2() results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.altf2 <- function(x, non.interactive=NULL, ...)
  {

if (is.null(non.interactive)) 
  {
    non.interactive <- FALSE
  }
  
nmods <- length(x$y.hat)

plot1g <- function(x)
  {
    
    col <- rich.colors(nmods, palette="temperature", rgb=FALSE, plot=FALSE) 

    names <- colnames(x$coeff.[[1]])
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
    
    width <-  480
    height <- 300
    

    m1 <- vector()
    m2 <- vector()
    
    for (j in 1:(length(names))) 
      {
        if (names(x$coeff.)[1]=="av. OLS")
          {
            m1.t <- matrix(x$coeff.[[1]][,j],ncol=1,nrow=length(x$y))
          }
        else
          {
            m1.t <- x$coeff.[[1]][,j]
          }
        for (ii in 1:nmods)
          {
            if (names(x$coeff.)[ii]=="av. OLS")
             {
               m1.t2 <- matrix(x$coeff.[[ii]][,j],ncol=1,nrow=length(x$y))
             }
            else
             {
               m1.t2 <- x$coeff.[[ii]][,j]
             }
            m1.t <- cbind(m1.t, m1.t2)
          }
       m1[j] <- min(m1.t,na.rm=TRUE)
       m2[j] <- max(m1.t,na.rm=TRUE)
     }


    for (j in 1:(length(names)))  
      {
        mypath <- file.path(getwd(), paste("altf2_coeff_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(5, 1, 2, 1))
        plot(index(x$y), rep(NA,length(index(x$y))), lty=2, type="l", col="white", ylim=c(m1[j],m2[j]), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:nmods)
          {
            if (names(x$coeff.)[ii]=="av. OLS")
              {
                plot(index(x$y),rep(x$coeff.[[ii]][,j],length(x$y)),col=col[ii],ylim=c(m1[j],m2[j]),
                     axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            else
              {
                plot(index(x$y),x$coeff.[[ii]][,j],col=col[ii],ylim=c(m1[j],m2[j]),
                     axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            par(new=TRUE)
          }
        axis(1, at=inc, labels=labs)
        legend('bottom', inset=c(0,-0.40), names(x$coeff.), lty=rep(1,ii), col=col[1:ii], ncol=4, cex=0.9) 
        dev.off()
      }


     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("altf2_coeff_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="altf2_coeff_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
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
    
    col <- rich.colors(nmods, palette="temperature", rgb=FALSE, plot=FALSE) 

    names <- colnames(x$coeff.[[1]])
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
    
    width <-  480
    height <- 300
    
    for (j in 1:(length(names)))  
      {
        mypath <- file.path(getwd(), paste("altf2_p_val_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(5, 1, 2, 1))
        plot(index(x$y), rep(NA,length(index(x$y))), lty=2, type="l", col="white", ylim=c(0,1), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:nmods)
          {
             if (names(x$p.val)[ii]=="av. OLS")
              {
                 plot(index(x$y),rep(x$p.val.[[ii]][,j],length(x$y)),col=col[ii],ylim=c(0,1),
                      axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            else
              {
                 if (! names(x$p.val)[ii]=="av. TVP")
                  {
                     plot(index(x$y),x$p.val.[[ii]][,j],col=col[ii],ylim=c(0,1),
                          axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
                  }
              }
            par(new=TRUE)
          }
        axis(1, at=inc, labels=labs)
        if ("av. TVP" %in% names(x$coeff.))
          {
            if (length(names(x$coeff.))==1)
              {
                legend('bottom', inset=c(0,-0.40), c(""), lty=rep(1,(ii-1)), col=col[1:(ii-1)], ncol=3, cex=0.9) 
              }
            else
              {
                legend('bottom', inset=c(0,-0.40), names(x$coeff.)[-length(names(x$coeff.))], lty=rep(1,(ii-1)), col=col[1:(ii-1)], ncol=3, cex=0.9) 
              }
          }
        else
          {
            legend('bottom', inset=c(0,-0.40), names(x$coeff.), lty=rep(1,ii), col=col[1:ii], ncol=3, cex=0.9) 
          }

        dev.off()
      }


     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("altf2_p_val_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="altf2_p_val_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
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

    names <- names(x$weights)
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
    
    width <-  480
    height <- 300
    
    for (j in 1:(length(names)))  
      {
        mypath <- file.path(getwd(), paste("altf2_weights_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(2, 1, 2, 1))
        plot(index(x$y), rep(NA,length(index(x$y))), lty=2, type="l", col="white", ylim=c(0,1), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:ncol(x$weights[[1]]))
          {
            if (names(x$weights)[j]=="av. OLS")
              {
                plot(index(x$y),rep(x$weights[[j]][,ii],length(x$y)),col=col[ii],ylim=c(0,1),
                     axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            else
              {
                plot(index(x$y),x$weights[[j]][,ii],col=col[ii],ylim=c(0,1),
                     axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            par(new=TRUE)
          }
        axis(1, at=inc, labels=labs)
        dev.off()
      }


     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("altf2_weights_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="altf2_weights_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
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
    
    col <- rich.colors(nmods, palette="temperature", rgb=FALSE, plot=FALSE) 

    names <- colnames(x$coeff.[[1]])
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
    
    width <-  480
    height <- 300
    
    for (j in 1:(length(names)))  
      {
        mypath <- file.path(getwd(), paste("altf2_rvi_", j, ".png", sep = ""))
        png(filename = mypath, height = height)
        par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(5, 1, 2, 1))
        plot(index(x$y), rep(NA,length(index(x$y))), lty=2, type="l", col="white", ylim=c(0,1), 
             axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
        par(new=TRUE)
        for (ii in 1:nmods)
          {
             if (names(x$rel.var.imp.)[ii]=="av. OLS")
              {
                 plot(index(x$y),rep(x$rel.var.imp.[[ii]][,j],length(x$y)),col=col[ii],ylim=c(0,1),
                      axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            else
              {
                 plot(index(x$y),x$rel.var.imp.[[ii]][,j],col=col[ii],ylim=c(0,1),
                      axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              }
            par(new=TRUE)
          }
        axis(1, at=inc, labels=labs)
        legend('bottom', inset=c(0,-0.40), names(x$coeff.), lty=rep(1,ii), col=col[1:ii], ncol=4, cex=0.9) 

        dev.off()
      }


     img <- list()
     for (i in 1:(length(names)))
      {
        mypath <- file.path(getwd(), paste("altf2_rvi_", i, ".png", sep = ""))
        img[[i]] <- readPNG(mypath)
      }

      png(filename="altf2_rvi_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
      par(mar=c(0,0,0,0))
      layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
      for(i in 1:(length(names))) 
        {
          plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
          rasterImage(img[[i]],0,0,1,1) 
        }
      dev.off()
  }


plot5g <- function(x)
  {
    
    col <- rich.colors(nmods, palette="temperature", rgb=FALSE, plot=FALSE) 

    names <- names(x$weights)
    
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
      {
        inc[i+1] <- floor(i * length(x$y)/7)
      }
    labs <- rownames(x$y)[inc]
    
    par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(7, 1, 2, 1))
    
    for (j in 1:(length(names)))  
      {
        if (names(x$exp.var.)[j]=="av. OLS")
          {
            plot(index(x$y),rep(as.vector(x$exp.var.[[j]]),length(x$y)),col=col[j],ylim=c(0,ncol(x$coeff.[[1]])),
                 axes=TRUE, xaxt='n', xlab='', ylab='', type="l", main='exp. var.')
            par(new=TRUE)
          }
        else
          {
            plot(index(x$y),as.vector(x$exp.var.[[j]]),col=col[j],ylim=c(0,ncol(x$coeff.[[1]])),
                 axes=TRUE, xaxt='n', xlab='', ylab='', type="l", main='exp. var.')
            par(new=TRUE)
          }
        axis(1, at=inc, labels=labs)
        legend('bottom', inset=c(0,-0.45), names(x$coeff.), lty=rep(1,j), col=col[1:j], ncol=4, cex=0.9) 
      }

  }


        if (non.interactive == FALSE) 
          {
            choices <- c("expected coefficients - separate plots (files in working directory)",
                         "p-values for t-tests - separate plots (files in working directory)",
                         "models' weights - separate plots (files in working directory)",
                         "relative variable importance (files in working directory)",
                         "expected number of variables (incl. constant)")
            pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
            switch(pick, plot1g(x), plot2g(x), plot3g(x), plot4g(x), plot5g(x))
          }
        else
          {
            plot1g(x)
            plot2g(x)
            plot3g(x)
            plot4g(x)
            plot5g(x)
          }
 
  }
  