"moderegbw" <- function(Y, X, method="CV-density", p.order=0, h1=NULL, h2=NULL, nstart = 4,
                         xinterval = quantile(X, probs=c(0.025, 0.975), names = FALSE),
                         df=5, ncomp=5, nboot=5){
  #########################################################################################
  # data structure
  #########################################################################################
  n = length(Y);
  h1ref = c(1.06*sd(X)*n^(-0.2));
  h2ref = c(1.06*sd(Y)*n^(-0.2));
  if(is.null(h1)) h1 = seq(h1ref*0.2, h1ref*1.5, length.out = 10); 
  if(is.null(h2)) h2 = h2ref;
  band.h = as.matrix(expand.grid(h1, h2));
  dd=density(X)
  dens=splinefun(dd$x, dd$y)
  fXhat=function(x) {
    resp = ifelse((x<=xinterval[2]&x>=xinterval[1]), 1, 0);
    resp[resp<0]=0; 
    resp
  }
  pX = fXhat(X);
  if(method=="CV-density"){
    if(p.order==0){
      CV = .Call("CVdens_LCfit", X, Y, pX, h1, h2, PACKAGE = "lpme")$CV; 
    }else{
      CV = .Call("CVdens_LLfit", X, Y, pX, h1, h2, PACKAGE = "lpme")$CV; 
    }
  }else if(method=="CV-mode"){
    dd = cbind(X, Y, pX); Xorder = order(dd[,1]);
    dd = dd[order(dd[,1]),];
    x = dd[,1];
    nx = length(x);
    nyx = nstart;
    x.num = rep(nyx, nx);
    yindx = c(0, cumsum(x.num));
    meshy = rep(0, yindx[nx+1]);
    quanprobs = seq(0.05, 0.95, length=nyx);
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
      meshy[indxl:indxr] = sort(yii)[round(length(indx)*quanprobs)];
    }
    mesh = cbind(rep(x,x.num), as.vector(meshy));
    if(p.order==0){
      CV = .Call("CVmode_LCfit", dd[,1], dd[,2], dd[,3], meshy, yindx, h1, h2, 200, 1e-4, PACKAGE = "lpme")$CV; 
    }else{
      CV = .Call("CVmode_LLfit", dd[,1], dd[,2], dd[,3], meshy, yindx, h1, h2, 200, 1e-4, PACKAGE = "lpme")$CV;
    }
  }else if(method=="bootstrap"){
    x = seq(min(xinterval), max(xinterval), length.out = 200); nx=length(x);
    wide = abs(x[2]-x[1]);
    y = seq(min(Y), max(Y), length.out = 1000); ny=length(y);
    ## Mesh points
    dd = data.frame(XX=X, YY=Y);
    dd = dd[order(dd[,1]),];
    nyx = nstart;
    x.num = rep(nyx, nx);
    yindx = c(0, cumsum(x.num));
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
      meshy[indxl:indxr] = sort(yii)[round(length(indx)*quanprobs)];
    }
    mesh = cbind(rep(x,x.num), as.vector((meshy)));
    # bootstrap settings
    basisX=splines::bs(X,df);
    fitall = flexmix::stepFlexmix(Y ~ basisX, nrep = 20, k = 1:ncomp)
    fit0 = flexmix::getModel(fitall, which = "AIC")
    param = flexmix::parameters(fit0); 
    coeff = param[-nrow(param),];
    sig = param[nrow(param),];
    mu = cbind(1,predict(basisX, newx=x))%*%coeff;
    wei = flexmix::prior(fit0);
    fy_x = matrix(0, nx, ny)
    for(i in 1:nx){
      tmp=rep(0, ny);
      for(j in 1:length(wei)){
        tmp = tmp + wei[j]*dnorm(y, mean=mu[i,j], sd=sig[j])
      }
      fy_x[i,] = tmp
    }
    dd=density(X)
    dens=splinefun(dd$x, dd$y)
    fXhat=function(x) {
      resp = dens(x); resp[resp<0]=0; resp
    }
    px = fXhat(x);
    # start bootstrap
    res = matrix(0, 3, nrow(band.h));
    for(j in 1:nrow(band.h)){
      ISE=rep(0, nboot);
      for(boo in 1:nboot){
        newY = (flexmix::rflexmix(fit0))$y[[1]][,1];
        fit = modereg(newY, X, mesh=mesh, bw=band.h[j,], p.order=p.order)$mode;
        dis = rep(0, nx);
        for(i in 1:nx){
          ymhat = fit[(yindx[i]+1):(yindx[i+1])];
          ym = y[.localMaxima(fy_x[i,])];
          tmp = rep(0, length(ym));
          for(ii in 1:length(ym)) tmp[ii] = min(abs(ymhat-ym[ii]));
          dis1 = max(tmp);
          tmp = rep(0, nyx);
          for(ii in 1:nyx) tmp[ii] = min(abs(ym-ymhat[ii]));
          dis2 = max(tmp);
          dis[i] = max(c(dis1, dis2))
        }
        ISE[boo]=sum( dis^2*px )*wide;
      }
      res[,j] = c(mean(ISE), band.h[j,1], band.h[j,2]);
    }
    CV = matrix(res[1,], length(h1), length(h2))
  }else{
    stop("This function only supports CV-density, CV-mode and bootstrap");
  }
  hhxy = as.vector(band.h[which.min(as.vector(CV)),]);
  output <- list(bw = hhxy,
                 CV = CV);
  class(output) <- c("moderegbw")
  output
}