
grid.tvp <- function (y,x,V,grid.lambda,W=NULL,kappa=NULL,parallel.grid=NULL,c=NULL)
{

### this is a wrapper of tvp()

### "foreach", "doParallel", "stats", and "xts" packages are required

### grid.lambda must be numeric vector 

#############################################

if (is.null(parallel.grid)) { parallel.grid <- FALSE }

#############################################

grid.lambda <- unique(grid.lambda)
grid.lambda <- sort(grid.lambda, decreasing=TRUE)

j.lambda <- NULL

if (parallel.grid == TRUE)
  {
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl)
    
    g.tvp <- foreach(j.lambda=grid.lambda, .packages=c("xts","fDMA")) %dopar%
              {
                tvp(y=y,x=x,V=V,lambda=j.lambda,W=W,kappa=kappa,c=c)
              }
          
    stopCluster(cl)
  }
else
  {
    g.tvp <- foreach(j.lambda=grid.lambda, .packages=c("xts","fDMA")) %do%
              {
                tvp(y=y,x=x,V=V,lambda=j.lambda,W=W,kappa=kappa,c=c)
              }
  }

err <- matrix(NA,nrow=length(grid.lambda),ncol=2)
rownames(err) <- grid.lambda
colnames(err) <- c("RMSE","MAE")

for (j in 1:length(grid.lambda))
  {
    x <- g.tvp[[j]]
    e <- (mean((as.vector(x$y)-as.vector(x$y.hat))^2,na.rm=TRUE))^0.5
    e <- round(e,digits=4)
    err[j,1] <- e
    e <- mean(abs(as.vector(x$y)-as.vector(x$y.hat)),na.rm=TRUE)
    e <- round(e,digits=4)
    err[j,2] <- e
  }

temp <- list(g.tvp,err)
names(temp) <- c("models","fq")
class(temp) <- "grid.tvp"
return(temp)

}
