
grid.roll.reg <- function (y,x=NULL,grid.window,parallel.grid=NULL,c=NULL)
{

### this is a wrapper of roll.reg()

### "foreach", "doParallel", "stats", and "xts" packages are required

### grid.window must be numeric vector 

#############################################

if (is.null(parallel.grid)) { parallel.grid <- FALSE }

#############################################

grid.window <- unique(grid.window)
grid.window <- sort(grid.window, decreasing=TRUE)

j.window <- NULL

if (parallel.grid == TRUE)
  {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    g.roll.reg <- foreach(j.window=grid.window, .packages=c("xts","fDMA")) %dopar%
              {
                roll.reg(y=y,x=x,window=j.window,c=c)
              }
          
    stopCluster(cl)
  }
else
  {
    g.roll.reg <- foreach(j.window=grid.window, .packages=c("xts","fDMA")) %do%
              {
                roll.reg(y=y,x=x,window=j.window,c=c)
              }
  }

err <- matrix(NA,nrow=length(grid.window),ncol=2)
rownames(err) <- grid.window
colnames(err) <- c("RMSE","MAE")

for (j in 1:length(grid.window))
  {
    x <- g.roll.reg[[j]]
    e <- (mean((as.vector(x$y)-as.vector(x$y.hat))^2,na.rm=TRUE))^0.5
    e <- round(e,digits=4)
    err[j,1] <- e
    e <- mean(abs(as.vector(x$y)-as.vector(x$y.hat)),na.rm=TRUE)
    e <- round(e,digits=4)
    err[j,2] <- e
  }

temp <- list(g.roll.reg,err)
names(temp) <- c("models","fq")
class(temp) <- "grid.roll.reg"
return(temp)

}
