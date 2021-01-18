
### it is suggested to close all graphical devices before plotting altf() results, i.e., to use the command graphics.off()

### requires "png" and "gplots" packages

plot.altf <- function(x, non.interactive=NULL, ...)
  {
if (is.null(non.interactive)) 
  {
    non.interactive <- FALSE
  }
  
nmods <- length(x$coeff.)

plot1g <- function(x)
  {
    
   mm1 <- vector()
   
   for (i in 1:nmods)
    {
      if (names(x$coeff.)[i] %in% c("rec. OLS","roll. OLS","TVP"))
        {
          mm1[i] <- 1
        }
      else
        {
          mm1[i] <- 0
        }
    }

   mm2 <- vector()
   n.c <- 0
   for (i in 1:nmods)
    {
      if (names(x$coeff.)[i] %in% c("TVP-AR(1)","TVP-AR(2)"))
        {
          mm2[i] <- 1
          n.c <- max(n.c,ncol(x$coeff.[[i]]))
        }
      else
        {
          mm2[i] <- 0
        }
    }

    col <- rich.colors((sum(mm1)+sum(mm2)), palette="temperature", rgb=FALSE, plot=FALSE) 
     
    plmods1 <- which(mm1==1)
    plmods2 <- which(mm2==1)
      
    names1 <- colnames(x$coeff.[[plmods1[1]]])
    names2 <- colnames(x$coeff.[[plmods2[1]]])
    if (n.c == 2)
      {
        names2 <- c("const","ar1")
      }
    if (n.c == 3)
      {
        names2 <- c("const","ar1","ar2")
      }
      
    inc <- vector()
    inc[1] <- 1
    for (i in 1:7)
     {
       inc[i+1] <- floor(i * length(x$y)/7)
     }
    labs <- rownames(x$y)[inc]
      
    width <-  480
    height <- 300

    min1 <- vector()
    max1 <- vector()
    
    min2 <- vector()
    max2 <- vector()

    if (sum(mm1)>0)
      {
        for (j in 1:length(names1)) 
          {
            m1.t <- x$coeff.[[plmods1[1]]][,j]

            for (ii in 1:sum(mm1))
              {
                m1.t2 <- x$coeff.[[plmods1[ii]]][,j]
                m1.t <- cbind(m1.t, m1.t2)
               }
            min1[j] <- min(m1.t,na.rm=TRUE)
            max1[j] <- max(m1.t,na.rm=TRUE)
          }
       }
       

    if (sum(mm2)>0)
      {
        for (j in 1:length(names2)) 
          {
            if (ncol(x$coeff.[[plmods2[1]]]) >= j) 
              { 
                m1.t <- x$coeff.[[plmods2[1]]][,j] 
              }
            else
              {
                m1.t <- x$coeff.[[plmods2[2]]][,j] 
              }

            for (ii in 1:sum(mm2))
              {
                if (ncol(x$coeff.[[plmods2[ii]]]) >= j)
                  {
                    m1.t2 <- x$coeff.[[plmods2[ii]]][,j]
                    m1.t <- cbind(m1.t, m1.t2)
                  }
               }
            min2[j] <- min(m1.t,na.rm=TRUE)
            max2[j] <- max(m1.t,na.rm=TRUE)
          }
      }
      
  if (sum(mm1)>0 && sum(mm2)==0)
    {
      names <- names1
      m1 <- min1
      m2 <- max1
      mtoplot <- sum(mm1)
      plmods <- plmods1
    }
  if (sum(mm1)==0 && sum(mm2)>0)
    {
      names <- names2
      m1 <- min2
      m2 <- max2
      mtoplot <- sum(mm2)
      plmods <- plmods2
    }
  if (sum(mm1)>0 && sum(mm2)>0)
    {
      names <- c(names1,names2[-1])
      m1 <- c(min1,min2[-1])
      m2 <- c(max1,max2[-1])
      m1[1] <- min(min1[1],min2[1],na.rm=TRUE)
      m2[1] <- max(max1[1],max2[1],na.rm=TRUE)
      mtoplot <- sum(mm1) + sum(mm2)
      plmods <- c(plmods1,plmods2)
    }
 
  if ( sum(mm1) > 0 || sum(mm2) > 0 ) 
      {
      
      for (j in 1:(length(names)))  
        { 
          mypath <- file.path(getwd(), paste("altf_coeff_", j, ".png", sep = ""))
          png(filename = mypath, height = height)
          par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(5, 1, 2, 1))
          plot(index(x$y), rep(NA,length(index(x$y))), lty=2, type="l", col="white", ylim=c(m1[j],m2[j]), 
               axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
          par(new=TRUE)
          for (ii in 1:mtoplot)
            {
              j.ind <- which(colnames(x$coeff.[[plmods[ii]]]) == names[j])
              if ( length(j.ind) > 0 )
                {
                  plot(index(x$y),x$coeff.[[plmods[ii]]][,j.ind],col=col[ii],ylim=c(m1[j],m2[j]),
                       axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
                  par(new=TRUE)
                }
            }
          axis(1, at=inc, labels=labs)
          if (sum(mm1)>0 && sum(mm2)>0)
            {
              if ( names[j] == c("const") )
                {
                  legend('bottom', inset=c(0,-0.40), names(x$coeff.)[plmods], lty=rep(1,ii), col=col[1:ii], ncol=5, cex=0.8) 
                }
              if ( ! names[j] == c("const") && names[j] %in% names1 )
                {
                  legend('bottom', inset=c(0,-0.40), names(x$coeff.)[plmods1], lty=rep(1,sum(mm1)), col=col[1:sum(mm1)], ncol=3, cex=0.8) 
                }
              if ( ! names[j] == c("const") && names[j] %in% names2 )
                {
                  legend('bottom', inset=c(0,-0.40), names(x$coeff.)[plmods2], lty=rep(1,sum(mm2)), col=col[sum(mm1)+1:sum(mm2)], ncol=2, cex=0.8) 
                }
            }
          else
            {
              legend('bottom', inset=c(0,-0.40), names(x$coeff.)[plmods], lty=rep(1,ii), col=col[1:ii], ncol=5, cex=0.8) 
            }
          dev.off()
        }
       
       

       img <- list()
       for (i in 1:(length(names)))
        {
          mypath <- file.path(getwd(), paste("altf_coeff_", i, ".png", sep = ""))
          img[[i]] <- readPNG(mypath)
        }

        png(filename="altf_coeff_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
        par(mar=c(0,0,0,0))
        layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
        for(i in 1:(length(names))) 
          {
            plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
            rasterImage(img[[i]],0,0,1,1) 
          }
        dev.off()   
      
      }
    
  }



plot2g <- function(x)
  {
    
   mm1 <- vector()
   
   for (i in 1:nmods)
    {
      if (names(x$coeff.)[i] %in% c("rec. OLS","roll. OLS"))
        {
          mm1[i] <- 1
        }
      else
        {
          mm1[i] <- 0
        }
    }

   if (sum(mm1)>0)
    {

      col <- rich.colors(sum(mm1), palette="temperature", rgb=FALSE, plot=FALSE) 
      
      plmods <- which(mm1==1)
      
      names <- colnames(x$coeff.[[plmods[1]]])
      
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
          mypath <- file.path(getwd(), paste("altf_p_val_", j, ".png", sep = ""))
          png(filename = mypath, height = height)
          par(xpd=TRUE, fig = c(0, 1, 0, 1), oma = c(2, 1, 1, 1), mar = c(5, 1, 2, 1))
          plot(index(x$y), rep(NA,length(index(x$y))), lty=2, type="l", col="white", ylim=c(0,1), 
               axes=TRUE, xaxt='n', xlab="", ylab="", main=names[j])
          par(new=TRUE)
          for (ii in 1:sum(mm1))
            {
              plot(index(x$y),x$p.val.[[plmods[ii]]][,j],col=col[ii],ylim=c(0,1),
                   axes=FALSE, xaxt='n', xlab='', ylab='', type="l", main='')
              par(new=TRUE)
            }
          axis(1, at=inc, labels=labs)
          legend('bottom', inset=c(0,-0.40), names(x$coeff.)[plmods], lty=rep(1,ii), col=col[1:ii], ncol=2, cex=0.9) 
          dev.off()
        }

       img <- list()
       for (i in 1:(length(names)))
        {
          mypath <- file.path(getwd(), paste("altf_p_val_", i, ".png", sep = ""))
          img[[i]] <- readPNG(mypath)
        }

        png(filename="altf_p_val_all.png", width = 2 * width, height = height * ceiling((length(names))/2))
        par(mar=c(0,0,0,0))
        layout(matrix(1:(2*ceiling((length(names))/2)), ncol=2, byrow=TRUE))
        for(i in 1:(length(names))) 
          {
            plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
            rasterImage(img[[i]],0,0,1,1) 
          }
        dev.off()   
    }
    
  }

        if (non.interactive == FALSE) 
          {
            choices <- c("expected coefficients - separate plots (files in working directory)",
                         "p-values for t-tests - separate plots (files in working directory)")
            pick <- menu(choices = paste(" ", choices), title = "\nMake a plot selection (or 0 to exit):")
            switch(pick, plot1g(x), plot2g(x))
          }
        else
          {
            plot1g(x)
            plot2g(x)
          }

 
  }
  