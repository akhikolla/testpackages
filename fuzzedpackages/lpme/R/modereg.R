"modereg" <- function(Y, W, bw, xgrid=NULL, sig=NULL, nstart=4, p.order=0, maxiter = 500, 
                      tol=.Machine$double.eps^0.25, mesh=NULL, PLOT=FALSE, ...){
  #########################################################################################
  # data structure
  #########################################################################################
  n = length(Y);
  sdW = sd(W); sdY = sd(Y);
  if(is.null(bw)) stop("please specify the bandwidth vector bw");
  if(is.null(mesh)){
    if(is.null(xgrid)){
      xgrid = seq(quantile(W, probs=0.025), quantile(W,probs=0.975), length.out = 100);
    }
    dd = cbind(W, Y);
    dd = dd[order(dd[,1]),];
    x = xgrid;
    nx = length(x);
    nyx = nstart;
    x.num = rep(nyx, nx);
    yindx = c(0, cumsum(x.num));
    #h_padding = (4/(2+4))^(1/(2+6))/n^(1/(2+6))*mean(apply(dd,2,sd));
    meshy = rep(0, yindx[nx+1]);
    quanprobs = seq(0.05, 0.95, length=nyx);
    if(nyx==1) quanprobs = 0.5;
    for(ii in 1:nx){
      indxl=yindx[ii]+1;
      indxr=yindx[ii+1];
      xdiff= dd[-1,1]-dd[-n,1];
      xdiff = xdiff[xdiff!=0];
      win0=quantile(abs(xdiff), 0.1);
      indx = which(dd[,1]>(x[ii]-win0*nyx)&dd[,1]<(x[ii]+win0*nyx));
      win = win0;
      while(length(indx)<10*nyx){
        win=win+win0;
        indx = which(dd[,1]>(x[ii]-win)&dd[,1]<(x[ii]+win));
      }
      yii = dd[indx, 2];
      #yii = yii[which(nn2(dd, cbind(x[ii], yii),k=1)$nn.dist<h_padding)];
      meshy[indxl:indxr] = sort(yii)[round(length(indx)*quanprobs)];
    }
    mesh = cbind(rep(x,x.num), as.vector((meshy)));
  }else{
    mesh = mesh[order(mesh[,1]),];
    x = unique(mesh[,1]); nx = length(x);
    x.num=as.vector(table(mesh[,1]));
    yindx = c(0, cumsum(x.num));
    meshy = mesh[,2];
  }
  
  if(is.null(sig)){
    if(p.order==0){
      fit = .Call("LCfitModeReg", x, meshy, yindx, W, Y, bw[1], bw[2],
                  maxiter, tol, PACKAGE = "lpme")$mode;
    }else{
      fit = .Call("LLfitModeReg", x, meshy, yindx, W, Y, bw[1], bw[2],
                  maxiter, tol, PACKAGE = "lpme")$mode;
    }
  }else{
    dt = 0.001;
    tt = seq(-1,1,dt);
    if(p.order==0){
      fit = .Call("LCfitModeRegLap2", x, meshy, yindx, W, Y, dt, tt,
                  sig, bw[1], bw[2], maxiter, tol, PACKAGE = "lpme")$mode
    }else{
      fit = .Call("LLfitModeRegLap2", x, meshy, yindx, W, Y, dt, tt,
                  sig, bw[1], bw[2], maxiter, tol, PACKAGE = "lpme")$mode
    }
  }
  
  ### pruning the fitted modes
  mesh_fitted = cbind(rep(x,x.num), as.vector((fit)));
  if(sum(is.na(fit))>0){
    na.indx = which(is.na(fit));
    other.indx = which(!is.na(fit));
    x.old = 1e100;
    first.na=TRUE;
    for(i in 1:length(na.indx)){
      x.value = mesh_fitted[na.indx[i],1];
      if(x.value!=x.old){
        ind.diff = TRUE;
      }else{
        ind.diff = FALSE;
      }
      x.old = x.value;
      y.x = mesh_fitted[which(mesh_fitted[,1]==x.value),2];
      na.indicator=is.na(y.x);
      na.num = sum(na.indicator);
      if(na.num>0){
        if(na.num==length(y.x)){
          if(first.na){
            message(paste("Warnings: The algorithm frails to converge at the following x values; please increase the bandwidths.", collapse = "\n"));
          }
          first.na=FALSE;
          fit[na.indx[i]] = 9999;
          if(ind.diff) cat(paste("x=", x.value, " ", sep=""));
        }else{
          fit[na.indx[i]] = (y.x[!na.indicator])[1];
        }
      }
    }
  }
  mesh_fitted = cbind(rep(x,x.num), as.vector((fit)));
  mesh_fitted_sd = cbind(rep(x,x.num), as.vector((fit))*sdW/sdY);
  
  ### clustering the fitted modes
  h_padding = (4/(2+4))^(1/(2+6))/n^(1/(2+6))*mean(apply(cbind(W, Y*sdW/sdY),2,sd));
  data.hclust = hclust(dist(mesh_fitted_sd), method="single");
  data.labels = cutree(data.hclust, h=h_padding);
  result = list()
  result$fitted = mesh_fitted
  result$labels = data.labels
  result$ordered = list()
  for(j in 1:max(data.labels)){
    w_tmp = which(data.labels==j)
    x_tmp = mesh_fitted[w_tmp,1]
    o_tmp = order(x_tmp)
    D_tmp = rbind(mesh_fitted[w_tmp,])[o_tmp,]
    result$ordered[[j]] = D_tmp
  }
  
  ## plot
  if(PLOT==TRUE){
    col0 = (1:max(data.labels))+1;
    plot(W,Y, pch=20, col="grey",...)
    points(mesh_fitted, col=col0[result$labels], pch=20)
  }
  
  #### Save to a list
  output <- list(xgrid = x, 
                 x.num = x.num, 
                 mode = fit,
                 mesh = mesh,
                 result = result);
  class(output) <- c("modereg")
  output
}
