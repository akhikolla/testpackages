#' @title ResautoDownscale
#'
#' @description This function is the iterative implementation of downscaling with autoencoder based residual network.
#'
#' @param r2 Stack for the covariates for downscaling   
#'
#' @param  fpredict0 Starting predictions obtained by GWR or other optimal learners such as XGBoost, residual network  
#'
#' @param  c.grid Coarsely resolved grid  
#'
#' @param  ss Sampling propotion for independent test  
#'
#' @param  nepoch Numder of epoch for residual network training (default: 30)
#'
#' @param  cores Number of CPU cores used for computing (default: 1)
#'
#' @param  thresh Stopping creation shreshold (default: 0.01) 
#'
#' @param  ntime Maximum number of iterations (default: 5) 
#'
#' @return List(performance metrics such as R2, RMSE, and downscaled images) 
#'
#' @export ResautoDownscale

ResautoDownscale=function(r2,fpredict0,c.grid,ss= 0.2,nepoch=30, cores= 1, thresh = 0.01,ntime=5){ 
  
  #r2,fpredict0,c.grid
  # ss= 0.2; cores= 5; thresh = 0.01; ntime=6
  metric_r2= keras::custom_metric("rsquared", function(y_true, y_pred) {
    SS_res =keras::k_sum(keras::k_square(y_true-y_pred ))
    SS_tot =keras::k_sum(keras::k_square(( y_true - keras::k_mean(y_true))))
    return ( 1 - SS_res/(SS_tot + keras::k_epsilon()))
  })
  lastpredict=fpredict0 
  validV=raster::values(fpredict0)
  validV=validV[!is.na(validV)]
  maxV=max(max(validV),mean(validV)+3*IQR(validV))
  minV=min(min(validV),mean(validV)-3*IQR(validV))
  minV=max(minV,0)
  diogRMSE= data.frame(tindex=0,stringsAsFactors = FALSE)  
  raster::beginCluster(cores)
  c.gridCELL=raster::rasterFromXYZ(cbind(raster::xyFromCell(c.grid, 
                  cell=seq(1,raster::ncell(c.grid)),spatial=FALSE),seq(1:raster::ncell(c.grid))))
  names(c.gridCELL)="cellNumbers_coarse"
  #resample coarse grids to fine grid
  c.grid_ds=raster::resample(c.grid, fpredict0,method="ngb")
  c.gridCELL_ds=raster::resample(c.gridCELL, fpredict0,method="ngb") 
  c.dat= raster::getValues(c.grid)
  c.datXY= raster::xyFromCell(c.grid, cell=seq(1,raster::ncell(c.grid)), 
                              spatial=FALSE)
  c.dat= cbind(seq(1,raster::ncell(c.grid)),c.datXY,c.dat)
  c.dat= as.data.frame(c.dat)
  names(c.dat)[1]= "cell"
  c.dat_ref= c.dat[which(complete.cases(c.dat)),]
  early_stopping = keras::callback_early_stopping(monitor ='loss', min_delta=0.000001)
  reduce=keras::callback_reduce_lr_on_plateau(patience=15)
  map1=lastpredict
  names(map1)="map1"
  for(i in c(1: (ntime))){ # i=1 
    message("Iteration ",i," ... ...") 
    r4= raster::stack(c.gridCELL_ds,map1)
    r4.dat=raster::getValues(r4)
    r4.dat= as.data.frame(r4.dat)
    r4.dat= r4.dat[which(complete.cases(r4.dat)),] #complete cases
    downFit=raster::aggregate(r4.dat$map1,list(group=r4.dat$cellNumbers_coarse),mean) # average within coarse grid
    xx= raster::merge(c.dat_ref, downFit, by.x = "cell", by.y = "group") 
    names(xx)= c("cell", "X", "Y", "C_grid", "F_grid")
    diogRMSE[i,"tindex"]=i 
    corr=cor(xx$C_grid,xx$F_grid)
    rsquare=rSquared(xx$C_grid,xx$C_grid-xx$F_grid)[1,1]
    rmse=rmse(xx$C_grid,xx$F_grid)
    varCalc= ((1)^2/(nrow(xx)*(nrow(xx)-1)))* sum((xx$C_grid -xx$F_grid)^2)
    seVar= sqrt(varCalc)
    tCrit= qt(1-(0.05/2), df=nrow(xx)-1)
    cI= seVar*tCrit
    upperCi= (mean((xx$C_grid -xx$F_grid)^2))+ cI
    lowerCi= (mean((xx$C_grid -xx$F_grid)^2))- cI
    diogRMSE[i,"lowerCi2"]= lowerCi
    diogRMSE[i,"upperCi2"]= upperCi
    diogRMSE[i,"mid"]=sqrt(mean((xx$C_grid -xx$F_grid)^2))
    
    diogRMSE[i,"resnetccr2"]=rsquare
    diogRMSE[i,"resnetccrmse"]=rmse
    diogRMSE[i,"resnetcccor"]=corr
    if(i==(ntime+1)){
      break 
    }
    if(i>3 && (diogRMSE[i-2,"resnetccr2"]-diogRMSE[i,"resnetccr2"])>0.02
       && (diogRMSE[i-2,"resnetccr2"]-diogRMSE[i-1,"resnetccr2"])>0.02){
      break;
    }
    xx$AF= -99999
    xx$AF[which(complete.cases(xx))]= xx$C_grid/xx$F_grid 
    xx[xx == -99999] = NA
    #make adjustment factor raster
    # fs<- paste("disseverOuts/",paste("iter_",zz,sep=""),sep="")
    AF.grid=raster::rasterFromXYZ(xx[,c(2,3,6)])
    #fine grid
    AF.grid_ds=raster::resample(AF.grid, fpredict0,method="ngb")
    #do the adjustment
    r5= raster::stack(map1,AF.grid_ds)
    f1 = function(x) (x[[1]]*x[[2]]) 
    upd.test = raster::clusterR(r5, fun=f1)
    names(upd.test)= names(c.grid_ds)
    ###DO the modelling
    #take a sample for modelling
    r1= raster::stack(upd.test,c.gridCELL_ds)  
    r3= raster::stack(r1,r2)
    vvfull=raster::values(r3)
    tlen=length(names(r3))
    scalev = raster::scale(vvfull[,c(3:tlen)]) 
    col_means = attr(scalev, "scaled:center") 
    col_stddevs = attr(scalev, "scaled:scale")
    vvfull[,c(3:tlen)]=raster::scale(vvfull[,c(3:tlen)], center = col_means, scale = col_stddevs)
    vv=na.omit(vvfull) # ss=0.2
    train_index=sample(c(1:nrow(vv)),size=ceiling(nrow(vv)*(1-ss)))
    test_index=setdiff(c(1:nrow(vv)),train_index)
    #sr=sampleRandom(r3,ceiling(ncell(r3[[1]])*ss)) # sample random grid cells
    sr= as.data.frame(vv[train_index,])
    x_train=as.matrix(sr[,c(3:tlen)])
    y_train=as.vector(as.matrix(sr[,1]))
    se= as.data.frame(vv[test_index,])
    x_test=as.matrix(se[,c(3:tlen)])
    y_test=as.vector(as.matrix(se[,1]))
    
    nfea=tlen-2;nout=1;nodes=c(32,16,8,4);mdropout=0.2;isres=TRUE;outtype=0;fact="linear"
    acts=rep("relu",length(nodes));fact="linear";reg=NULL;batchnorm=TRUE
    autoresmodel=AutoEncoderModel(nfea,nout,nodes,acts,mdropout,reg,batchnorm,isres,outtype,fact=fact)
    autoresmodel %>% keras::compile(
      loss = "mean_squared_error",
      optimizer = keras::optimizer_rmsprop(),
      metrics = c("mean_squared_error",metric_r2)
    )
    history <- autoresmodel %>% keras::fit(
      x_train, y_train, 
      epochs = nepoch, batch_size = 2560, callbacks=list(early_stopping,reduce),
      validation_split = 0.2,verbose=0
    )
    pred = autoresmodel %>% predict(x_test,batch_size = 2560)
    rquared=rSquared(y_test,y_test-pred)[1,1]
    rmse=rmse(y_test,pred)
    na_index=which(!complete.cases(vvfull[,c(3:tlen)]))
    predfull = autoresmodel %>% predict(vvfull[,c(3:tlen)])
    predfull[which(predfull<minV)]=minV
    predfull[which(predfull>maxV)]=maxV
    predfull[na_index]=NA 
    findex=length(history[["metrics"]]$rsquared) 
    diogRMSE[i,"resnetr2"]=history[["metrics"]]$rsquared[[findex]]
    diogRMSE[i,"resnetrmse"]=history[["metrics"]]$mean_squared_error[[findex]]
    diogRMSE[i,"resnetvalr2"]=history[["metrics"]]$val_rsquared[[findex]]
    diogRMSE[i,"resnetvalrmse"]=history[["metrics"]]$val_mean_squared_error[[findex]]
    diogRMSE[i,"resnettestr2"]=rquared 
    diogRMSE[i,"resnettestrmse"]=rmse 
    map2=fpredict0
    raster::values(map2)=predfull 
    names(map2)= "map2" 
#    map2=mask(map2,fpredict0)
    message("Modeling RMSE = ",round(rmse,3),";r2=",rquared)
    if (i >= 3) {
      FF= mean(abs(diogRMSE[i-2,"mid"]-diogRMSE[i-1,"mid"])+ 
                 abs(diogRMSE[i-1,"mid"]-diogRMSE[i,"mid"]) + abs(diogRMSE[i-2,"mid"]-diogRMSE[i,"mid"]))
      if (FF <= thresh) {break}
    } 
    map1= map2
    names(map1)= "map1"
    if(i==1){
      maxr2=rquared
      minrmse=rmse 
      maxmap=map1
    }else if(i>1 && maxr2<rquared){
      maxr2=rquared
      maxmap=map1
      minrmse=rmse 
    }
    rm(list=c("predfull","map2","autoresmodel","vvfull","vv","r1","r3","scalev",
              "x_train","y_train","x_test","y_test","sr","se","r4.dat","r4",
              "train_index","test_index"))
    gc()
  }
  #names(diogRMSE)=c("tindex","lowerCi","upperCi","mid","resnetccr2","resnetccrmse","resnetcccor",
  #                  "resnetr2","resnetrmse","resnetvalr2","resnetvalrmse","resnettestr2","resnettestrmse")
  raster::endCluster() 
  return(list(diogRMSE=diogRMSE,raster=maxmap,r2=maxr2,rmse=minrmse))
}